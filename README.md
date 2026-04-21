# Higgs Trilinear Self-Coupling (λ₃) via ZH Process

## Part 1: Physics Background

The Higgs trilinear self-coupling λ₃ (the HHH vertex) cannot be directly measured in single-Higgs production processes (like ZH) at tree level. However, λ₃ enters these processes **at one loop** through virtual EW corrections involving the H→HH splitting.

By calculating the **ratio of the NLO EW correction to the LO cross section** — parameterized as the **C₁** coefficient — we can extract the specific sensitivity of the process to the true value of λ₃.
For example, for the ZH process at 13 TeV across the inclusive phase space, this factor is roughly **C₁ ≈ 1.19%**.

The full BSM cross section for arbitrary κ_λ (where $\kappa_\lambda = \lambda_3 / \lambda_3^{SM}$) is reconstructed identically on an event-by-event level using analytic rescaling:

$$ \sigma(\kappa_\lambda) = \sigma_{LO} \times Z_H^{BSM} \times \left[ 1 + \kappa_\lambda \cdot C_1 + \delta_{ZH} + \dots \right] $$

where:
- $Z_H^{BSM} = \frac{1}{1 - (\kappa_\lambda^2 - 1) \delta_{ZH}}$
- $\delta_{ZH} \approx -1.536 \times 10^{-3}$
- $C_1 = \frac{\Delta \sigma_{EW}}{\sigma_{LO}}$

The pipeline automates this by generating LO QCD events, dynamically reweighting them to compute the exact NLO EW $O(\lambda_3)$ corrections for the corresponding physical phase space point, and then analytically applying the model formula above to compute final inclusive weights corresponding to arbitrary values of $\kappa_\lambda$.

---

## Quick Start

Three commands cover the normal workflow:

```bash
./setup_env.sh
./run_pipeline.sh --nevents 500000 --ecm 13600
./analyze.sh --kappa 0,1,2,5,10,-2,-5,-10
```

- `setup_env.sh` checks for Apptainer, installs LHAPDF 6.2.1 and MG5_aMC_v2.5.5 if missing, configures the UFO model, and prepares PDF libraries.
- `run_pipeline.sh` generates the MG5 processes, builds the OLP reweighter, produces LO events, and reweights them to attach EW loop corrections. Results are copied to `output/`.
- `analyze.sh` converts the LHE files to ROOT, computes BSM weights for each requested $\kappa_\lambda$, and runs the full plotting suite. All plots land in `output/plots/`.

---

## Prerequisites

- Linux environment with `bash`
- `apptainer` (or `singularity`)
- A local container image. The pipeline defaults to `trilinear-boost.sif` (build from `docker://lordcristo/trilinear-boost` if needed):
  ```bash
  apptainer build trilinear-boost.sif docker://lordcristo/trilinear-boost
  ```
- Python 3 with packages from `requirements.txt`:
  ```bash
  python3 -m pip install -r requirements.txt
  ```

**Required Python packages:**
- `uproot>=5.0`
- `numpy`
- `matplotlib`
- `awkward>=2.0`

---

## Script Options

### `run_pipeline.sh`

| Option | Default | Description |
|---|---:|---|
| `--nevents N` | `500000` | Number of LO events to generate |
| `--ecm GeV` | `13600` | Proton-proton center-of-mass energy |
| `--skip-generation` | — | Reuse an existing `hz_MC/` directory (skip MG5 process generation) |

### `analyze.sh`

| Option | Default | Description |
|---|---:|---|
| `--kappa LIST` | `0,1,2,5,10,-2,-5,-10` | Comma-separated $\kappa_\lambda$ values to process |
| `--ebeam GeV` | `6800` | Beam energy per proton passed to the LHE→ROOT conversion |
| `--input-dir DIR` | `output/` | Directory containing `events.lhe` and `events_rwgt.lhe` |

---

## Output Layout

All generated artifacts live under `output/`:

```text
output/
├── events.lhe              # LO events from MG5
├── events_rwgt.lhe         # Reweighted events with NLO EW loop corrections
├── events_lo.root          # LO events in ROOT format
├── events_rwgt.root        # Reweighted events in ROOT format
├── events_lo_l3corr.root   # LO events corrected with SM NLO weights
├── events_l3corr_*.root    # BSM-weighted ROOT files for each kappa
└── plots/
    ├── weight_ratio.png
    ├── C1_vs_pt.png
    ├── compare_and_C1.png
    ├── compare_scatter.png
    ├── compare_hpt_hist.png
    ├── compare_zpt_hist.png
    └── events_l3corr_*_vs_lo.png
```

`plots/` contains all PNG output from the validation and plotting scripts.

---

## Manual Workflow / Debug Reference

For step-by-step debugging without the wrapper scripts, see [`manual_stage3_debug.md`](manual_stage3_debug.md).

The manual flow is:

1. **Environment** — source `lhapdf_env.sh` and enter the Apptainer container.
2. **MG5 LO generation** — run `./bin/mg5_aMC` with `cards/proc_mc_hz` to create `hz_MC`.
3. **Virtual extraction** — run `trilinear-RW/gevirt.sh hz_MC` to produce `check_olp.inc` and `proc_ml`.
4. **ME generation** — generate `hz_ME` with the loop-diagram filter (`vvh-loop_diagram_generation.py`) and `hhh-model-new`.
5. **Build OLP** — compile `check_OLP` in `hz_ME/SubProcesses/`.
6. **LO events** — launch MG5 for LO event generation (`noshowerLO`).
7. **Reweighting** — unzip `events.lhe.gz` and run `./check_OLP` to produce `events_rwgt.lhe`.
8. **Analysis** — run `lhe_to_root.py`, `add_l3_weight.py`, and the plotting scripts.

The wrapper scripts (`run_pipeline.sh` and `analyze.sh`) apply several automatic fixes that the manual guide documents in detail, including Boost-header patches for `pdf_lhapdf6.cc` and the `pdlabel='nn23nlo'` fallback in `check_OLP.f`.
