# VH Workflow Extension Design

Date: 2026-05-16

## Goal

Extend the current `VHtrilinear` repository from a single-process `ZH` pipeline into a shared `VH` workflow that supports:

- `ZH` as an isolated end-to-end process
- `WH` as an isolated end-to-end process
- `WH` implemented internally as two subchannels, `W+H` and `W-H`
- only combined `WH` products exposed as public outputs

The user-facing interface should be generic, but the process execution and artifacts must remain isolated from beginning to end.

## Scope

This design covers:

- process selection in setup, generation, and analysis wrappers
- process-specific card layout
- process-specific MG5 and reweighting work directories
- process-specific output layout under `output/`
- internal `WH` subchannel execution and merge
- process-aware analysis and plotting

This design does not cover:

- adding `VBF`, `tHj`, or `ttH`
- changing the upstream `trilinear-RW` physics kernel
- retuning process-specific constants beyond what is needed to support `WH`

## Requirements

1. The external interface must support process selection via a flag such as `--process zh|wh`.
2. `ZH` and `WH` must not share generated process directories, intermediate reweighting products, or final output trees.
3. `WH` must internally run `W+H` and `W-H` separately.
4. Public `WH` outputs must be combined products only.
5. Combined `WH` ROOT outputs must preserve a `subchannel` branch so `W+H` and `W-H` contributions can be audited later.
6. Weighted yields must remain physically correct; the design must not assume raw event counts are proportional to cross section.
7. Validation and plot labeling must be process-aware so `ZH` labels do not appear in `WH` output.

## Recommended Approach

Use a single shared `VH` driver with a small process-spec layer.

The process spec defines, for each public process:

- process key: `zh` or `wh`
- displayed label: `ZH` or `WH`
- public output root
- process cards
- internal subchannels
- MG5 workdir names
- plotting labels
- merge behavior

This keeps the orchestration code shared while making all process-specific state explicit and isolated.

## Interface

### Wrapper scripts

The wrappers become:

- `./run_pipeline.sh --process zh`
- `./run_pipeline.sh --process wh`
- `./analyze.sh --process zh`
- `./analyze.sh --process wh`

Optional flags like `--nevents`, `--ecm`, and analysis options remain available.

### Process semantics

- `zh` means one public process with one internal generation path.
- `wh` means one public process with two hidden internal subchannels:
  - `wh_plus` for `p p > h w+`
  - `wh_minus` for `p p > h w-`

Users do not run `wh_plus` or `wh_minus` directly through the public interface.

## Repository Structure

### Cards

Move process inputs under process-specific paths:

- `cards/zh/proc_card.dat`
- `cards/zh/run_card.dat`
- `cards/zh/param_card.dat`
- `cards/wh/proc_card_wp.dat`
- `cards/wh/proc_card_wm.dat`
- `cards/wh/run_card.dat`
- `cards/wh/param_card.dat`

The `WH` run and parameter cards can initially be shared between `W+H` and `W-H` unless MG5 or physics requirements force divergence.

### Outputs

Public outputs are separated:

- `output/zh/`
- `output/wh/`

Examples:

- `output/zh/events.lhe`
- `output/zh/events_rwgt.lhe`
- `output/zh/events_lo.root`
- `output/zh/events_rwgt.root`
- `output/zh/plots/...`

- `output/wh/events_lo.root`
- `output/wh/events_rwgt.root`
- `output/wh/events_l3corr_*.root`
- `output/wh/plots/...`

Internal `WH` subchannel artifacts should live in process-private temporary locations and should not be treated as public deliverables.

### MG5 workdirs

Generated workdirs must be process-qualified.

Examples:

- `zh_MC`, `zh_ME`
- `whp_MC`, `whp_ME`
- `whm_MC`, `whm_ME`

This avoids cross-talk between `ZH` and `WH`, and also between `W+H` and `W-H`.

## Execution Flow

### ZH

`ZH` keeps the current physics flow:

1. generate LO process
2. run `gevirt.sh`
3. generate EW virtual process
4. build `check_OLP`
5. generate LO events
6. run reweighting
7. convert LHE to ROOT
8. compute `kappa_lambda` weights
9. produce plots

The main change is process qualification of paths and labels.

### WH

`WH` runs the same flow twice:

1. run full pipeline for `wh_plus`
2. run full pipeline for `wh_minus`
3. convert both subchannels to ROOT
4. merge the ROOT outputs into combined `WH` products
5. run weight construction and plotting on the merged ROOT products

The merge point is after LHE-to-ROOT conversion, not at raw LHE level.

## Merge Strategy

### Why merge at ROOT level

Merge after `lhe_to_root.py` because this is the simplest reliable boundary:

- the schema is explicit
- event counts and summed weights are easy to validate
- no custom LHE concatenation logic is needed
- downstream analysis already consumes ROOT

### Merge semantics

The combined `WH` ROOT products are formed by concatenating the `W+H` and `W-H` event tables.

The merge must preserve:

- all existing kinematic branches
- the `weight` branch
- `event_id`
- a new `subchannel` branch

`subchannel` values should be stable string or integer-coded labels that distinguish:

- `wp`
- `wm`

### Event normalization

If `W+H` and `W-H` are generated with equal requested `nevents`, the row counts will be roughly 50:50, but the weighted yields will still be physical because each event carries the MG5 normalization in its weight.

The design must treat weighted sums, not raw row counts, as the physics yield.

No artificial rescaling to force a target `W+H`:`W-H` event-count ratio is required.

## Analysis Design

### LHE to ROOT

The existing converter already recognizes either `Z` or `W` as the vector boson. It should be formalized as process-aware rather than implicitly generic.

Expected changes:

- allow process metadata or subchannel metadata to be attached during conversion
- write `subchannel` for merged `WH` products
- keep branch naming generic as `v_*`, not `z_*`

### Weight construction

The current `add_l3_weight.py` logic stays shared. It operates on ROOT products, so after `WH` merge it can run once on the combined `WH` LO and reweighted ROOT files.

The process config should supply process-specific labels and any process-specific physics constants if `WH` needs values different from `ZH`.

### Plotting

Plotting scripts need a process label parameter so hard-coded `ZH` titles are removed from shared code.

At minimum:

- process name in figure text must be configurable
- axis labels that say `Z` specifically must become `V` or process-aware
- output locations must use `output/<process>/plots/`

## Validation

### ZH validation

Retain the current checks, but report them within the `output/zh/` context.

### WH validation

Validate both subchannels before merge:

- expected files exist for `wh_plus`
- expected files exist for `wh_minus`
- ROOT schemas match
- subchannel event counts are nonzero
- summed weights are finite

Validate merged `WH` outputs after merge:

- merged event count equals the sum of both subchannels
- merged weighted sums equal the sum of subchannel weighted sums
- `subchannel` branch exists
- both `wp` and `wm` appear

### Reporting

For `WH`, logs should report both:

- raw event counts
- weighted sums

This prevents confusion between sampling statistics and physics yield.

## Failure Handling

1. If one `WH` subchannel fails, the full `WH` run fails.
2. Merge must fail on schema mismatch.
3. Analysis must fail if process-specific inputs are missing.
4. Process-specific temporary paths must not overwrite another process run.
5. Plotting must fail fast if process metadata is inconsistent.

## Implementation Notes

The existing code is close enough that this is an adaptation, not a rewrite.

The likely refactor centers are:

- `run_pipeline.sh`: replace `hz_*` hard-coding with process/subchannel specs
- `analyze.sh` and `scripts/analyze.py`: route to process-specific inputs and outputs
- `scripts/lhe_to_root.py`: preserve process-aware metadata and support `subchannel`
- plotting scripts: remove hard-coded `ZH` labels
- cards: split into `zh` and `wh` process directories

## Tradeoffs

### Chosen tradeoff

Keep one shared driver with explicit process specs.

Benefits:

- less duplicated orchestration code
- easier to extend safely later
- keeps `ZH` and `WH` isolated without creating two drifting code paths

Cost:

- requires a careful first refactor of path handling and labels

### Rejected tradeoff

Copy the `ZH` flow into a separate `WH` implementation.

Reason:

- fastest short-term path
- highest long-term maintenance cost
- almost guaranteed drift in fixes and behavior

## Success Criteria

The design is successful when:

1. `./run_pipeline.sh --process zh` still works.
2. `./run_pipeline.sh --process wh` produces a combined `WH` result from internal `W+H` and `W-H` runs.
3. `./analyze.sh --process zh` and `./analyze.sh --process wh` write outputs under separate trees.
4. `output/wh/` contains only combined public artifacts.
5. Combined `WH` ROOT outputs contain a `subchannel` branch for auditing.
6. Weighted yields remain correct for merged `WH` products.
7. Shared plotting and analysis code no longer leaks `ZH`-specific labels into `WH`.
