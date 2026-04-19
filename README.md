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

## Part 2: Operating Instructions

This section outlines the steps to build and run the entire computational pipeline, from MadGraph event generation, matrix manipulation, to final ROOT-based statistical extraction and plotting.

### 1. Environment Setup

The pipeline requires MG5_aMC **v2_5_5**, LHAPDF 6.2.1, and an Apptainer container supporting GCC ≥ 4.6.0. 

Always align the environment pathings to standard configurations before starting:
```bash
# From the project root directory
source lhapdf_env.sh
```

### 2. MadGraph Process Generation (LO QCD)

Generate the LO background utilizing the custom `hhh-model` UFO.

```bash
cd MG5_aMC_v2_5_5
./bin/mg5_aMC
```
In the MadGraph prompt:
```madgraph
import model hhh-model
generate p p > h z [LOonly= QCD]
output hz_MC
quit
```

### 3. Virtual EW Subprocesses Extraction

Use the provided shell script to extract the required EW matrix elements.
```bash
cp ../trilinear-RW/gevirt.sh ./
./gevirt.sh hz_MC
```
This produces the Fortran include `check_olp.inc` and process target `proc_ml`. 

Next, generate the EW matrix elements loop structures. *(Ensure **COLLIER** is thoroughly disabled in your `input/mg5_configuration.txt`, by setting `collier = None`)*:

```bash
# Setup the filter script
cp ../trilinear-RW/vvh-loop_diagram_generation.py madgraph/loop/loop_diagram_generation.py

# In MadGraph prompt:
import model hhh-model
# (copy/paste the contents of proc_ml here line by line)
output hz_ME
quit
```

### 4. Build the Reweighting Executable

You will need to manually compile the OLP (One Loop Provider) evaluation code.
```bash
cd hz_ME/SubProcesses/

# Copy necessary files from across the ecosystem
cp ../../../trilinear-RW/makefile .
cp ../../../trilinear-RW/check_OLP.f .
cp ../../check_olp.inc .
cp P0_uux_hz/pmass.inc .
cp P0_uux_hz/nsqso_born.inc .
cp P0_uux_hz/nsquaredSO.inc .
cp ../../hz_MC/SubProcesses/c_weight.inc .
cp ../../hz_MC/SubProcesses/P0_uux_hz/nexternal.inc .

# Setup PDF libraries
cd ../lib/
cp ../../lhapdf/lib/libLHAPDF.a .
cp -r ../../lhapdf/share/LHAPDF PDFsets  # Link or copy LHAPDF sets

# Build the executable
cd ../SubProcesses/
make OLP_static
make check_OLP
```
*(This produces the critical executable `./check_OLP` taking raw events matching topology and assigning complex matrix components)*

### 5. Event Generation and Reweighting

1. **LO Events:** In `hz_MC/Cards/run_card.dat`, set `True = store_rwgt_info`. Launch the `hz_MC` process ensuring `fixedorder=OFF`, `shower=OFF`, `reweight=OFF`, `order=LO`. *(Generates `events.lhe.gz`)*
2. **Reweighting:** Unzip the generated LHE file and feed it to the OLP executable to attach the EW loop weights.
```bash
gunzip -c hz_MC/Events/run_01_LO/events.lhe.gz > hz_ME/SubProcesses/events.lhe
cd hz_ME/SubProcesses/
./check_OLP
```
*(Outputs `events_rwgt.lhe` which includes explicit NLO loop weights directly mapped).*

### 6. LHE to ROOT Conversion and BSM Weights computation

Leverage the internal analysis scripts to parse the text LHE formats into performant ROOT trees.

```bash
# 1. Convert LHE to flat ROOT Branches
python3 scripts/lhe_to_root.py events.lhe scripts/events_lo.root
python3 scripts/lhe_to_root.py events_rwgt.lhe scripts/events_rwgt.root

# 2. Append generic kappa_lambda event weights
# Compute C1 implicitly and output model weights. Example SM limit (kappa_lambda = 1.0)
python3 scripts/add_l3_weight.py scripts/events_lo.root scripts/events_rwgt.root --l3 1.0

# Generate templates for scale anomalous variants (e.g. kappa_lambda = -5.0)
python3 scripts/add_l3_weight.py scripts/events_lo.root scripts/events_rwgt.root --l3 -5.0
```
*(Outputs suffix specific roots evaluating multiple cross section states simultaneously `events_rwgt_l3corr_-5p0.root`... )*

### 7. Validation & Plotting Suite

Several scripts are embedded mapped directly to standard observables (`pt`, `eta`, `cos_theta_star`) enabling pipeline auditing:

```bash
# Check C1 inclusive ratio map between generated events
python3 scripts/plot_weight_ratio.py scripts/events_lo.root scripts/events_rwgt.root

# Check C1 systematic differential variation across Higgs Transverse Momentum
python3 scripts/plot_C1_vs_pt.py scripts/events_lo.root scripts/events_rwgt.root

# Draw 2-panel validation maps revealing relative shape divergences 
python3 scripts/compare_and_C1.py scripts/events_lo.root scripts/events_rwgt.root

# Visually construct multiscale BSM samples vs Nominal distributions overlay profiles
python3 scripts/plot_kappa3.py scripts/events_rwgt_l3corr.root scripts/events_lo.root --feature h_pt
```