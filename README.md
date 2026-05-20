# Higgs Trilinear Self-Coupling (λ₃) via VH Process

## Part 1: Physics Background

The Higgs trilinear self-coupling λ₃ (the HHH vertex) cannot be directly measured in single-Higgs production processes (like ZH or WH) at tree level. However, λ₃ enters these processes **at one loop** through virtual EW corrections involving the H→HH splitting.

By calculating the **ratio of the NLO EW correction to the LO cross section** — parameterized as the **C₁** coefficient — we can extract the specific sensitivity of the process to the true value of λ₃.
For example, for the VH process at 13 TeV across the inclusive phase space, this factor is roughly **C₁ ≈ 1.19%**.

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
./run_pipeline.sh --process zh --nevents 500000 --ecm 13600
./analyze.sh --process zh --kappa 0,1,2,5,10,-2,-5,-10

./run_pipeline.sh --process wh --nevents 500000 --ecm 13600
./analyze.sh --process wh --kappa 0,1,2,5,10,-2,-5,-10
```

- `setup_env.sh` checks for Apptainer, installs LHAPDF 6.2.1 and MG5_aMC_v2.5.5 if missing, configures the UFO model, and prepares PDF libraries.
- `run_pipeline.sh` generates the MG5 processes, builds the OLP reweighter, produces LO events, and reweights them to attach EW loop corrections. Results are copied to `output/<process>/`.
- `analyze.sh` converts the LHE files to ROOT, computes BSM weights for each requested $\kappa_\lambda$, and runs the full plotting suite. All plots land in `output/<process>/plots/`.

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
| `--process` | `zh` | Public process to generate (`zh` or `wh`) |
| `--nevents N` | `500000` | Number of LO events to generate |
| `--ecm GeV` | `13600` | Proton-proton center-of-mass energy |
| `--skip-generation` | — | Reuse an existing MC directory (skip MG5 process generation) |

### `analyze.sh`

| Option | Default | Description |
|---|---:|---|
| `--process` | `zh` | Public process to analyze (`zh` or `wh`) |
| `--kappa LIST` | `0,1,2,5,10,-2,-5,-10` | Comma-separated $\kappa_\lambda$ values to process |
| `--ebeam GeV` | `6800` | Beam energy per proton passed to the LHE→ROOT conversion |

---

## Machine Learning Pipeline (C₁ Regressor)

To apply the C₁ correction directly to reconstructed events (e.g., NanoAOD) without requiring OLP reweighting, we train an XGBoost regressor to predict C₁ purely from event kinematics.

### 1. Model Training
Train the XGBoost model using the generated LO and reweighted ROOT files:
```bash
python3 scripts/train_c1_regressor.py --process zh
python3 scripts/train_c1_regressor.py --process wh
```
This generates the trained model (`output/<process>/c1_regressor/c1_regressor.json`), an ONNX export (`c1_regressor.onnx`), and training validation profile plots.

### 2. NanoAOD Prediction & Validation
Apply the trained regressor to NanoAOD files using their `LHEPart` branches:
```bash
python3 scripts/predict_c1_nano.py \
    --process zh \
    --input nanoAOD_temp/*.root \
    --validate

python3 scripts/predict_c1_nano.py \
    --process wh \
    --input nanoAOD_temp/*.root \
    --validate
```
The `--validate` flag performs a closure test, generating kinematic profile plots comparing the NanoAOD predictions back against the ground-truth LO sample in `output/<process>/plots/nano_validation/`.

---

## Output Layout

All generated artifacts live under `output/`:

```text
output/
├── zh/
│   ├── events.lhe
│   ├── events_rwgt.lhe
│   ├── events_lo.root
│   ├── events_rwgt.root
│   ├── events_l3corr_*.root
│   ├── c1_regressor/
│   │   ├── c1_regressor.json
│   │   ├── c1_regressor.onnx
│   │   └── *.png
│   ├── nano_c1_predictions.root
│   └── plots/
└── wh/
    ├── events_lo.root
    ├── events_rwgt.root
    ├── events_l3corr_*.root
    ├── c1_regressor/
    │   ├── c1_regressor.json
    │   ├── c1_regressor.onnx
    │   └── *.png
    ├── nano_c1_predictions.root
    └── plots/
```

`plots/` contains all PNG output from the validation and plotting scripts.

---

## Manual Workflow / Debug Reference

For step-by-step debugging without the wrapper scripts, see [`manual_stage3_debug.md`](manual_stage3_debug.md).

The manual flow is:

1. **Environment** — source `lhapdf_env.sh` and enter the Apptainer container.
2. **MG5 LO generation** — run `./bin/mg5_aMC` with `cards/zh/proc_card.dat` (or `cards/wh/proc_card_wp.dat` / `cards/wh/proc_card_wm.dat`) to create the MC process directory.
3. **Virtual extraction** — run `trilinear-RW/gevirt.sh <mc_dir>` to produce `check_olp.inc` and `proc_ml`.
4. **ME generation** — generate the ME process directory with the loop-diagram filter (`vvh-loop_diagram_generation.py`) and `hhh-model-new`.
5. **Build OLP** — compile `check_OLP` in the ME `SubProcesses/` directory.
6. **LO events** — launch MG5 for LO event generation (`noshowerLO`).
7. **Reweighting** — unzip `events.lhe.gz` and run `./check_OLP` to produce `events_rwgt.lhe`.
8. **Analysis** — run `lhe_to_root.py`, `add_l3_weight.py`, and the plotting scripts.

The wrapper scripts (`run_pipeline.sh` and `analyze.sh`) apply several automatic fixes that the manual guide documents in detail, including Boost-header patches for `pdf_lhapdf6.cc` and the `pdlabel='nn23nlo'` fallback in `check_OLP.f`.
