#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# VHtrilinear Stage 4: Analysis wrapper
# ---------------------------------------------------------------------------

SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPTDIR"

# Defaults
KAPPA="0,1,2,5,10,-2,-5,-10"
EBEAM=6800
INPUTDIR="output"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --kappa)
            KAPPA="$2"
            shift 2
            ;;
        --ebeam)
            EBEAM="$2"
            shift 2
            ;;
        --input-dir)
            INPUTDIR="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--kappa 0,1,2,...] [--ebeam 6800] [--input-dir output/]"
            exit 1
            ;;
    esac
done

# 1. Ensure output/plots exists
mkdir -p output/plots

# 2. Check LHE inputs
EVENTS_LO_LHE="${INPUTDIR}/events.lhe"
EVENTS_RWGT_LHE="${INPUTDIR}/events_rwgt.lhe"

if [[ ! -f "$EVENTS_LO_LHE" ]]; then
    echo "ERROR: $EVENTS_LO_LHE not found"
    exit 1
fi
if [[ ! -f "$EVENTS_RWGT_LHE" ]]; then
    echo "ERROR: $EVENTS_RWGT_LHE not found"
    exit 1
fi

# 3. Convert LHE -> ROOT
LO_ROOT="output/events_lo.root"
RWGT_ROOT="output/events_rwgt.root"

echo "[analyze] Converting LHE to ROOT (ebeam=${EBEAM} GeV)..."
python3 scripts/lhe_to_root.py "$EVENTS_LO_LHE" "$LO_ROOT" --ebeam "$EBEAM"
python3 scripts/lhe_to_root.py "$EVENTS_RWGT_LHE" "$RWGT_ROOT" --ebeam "$EBEAM"

# 4. Loop over kappa values and compute l3 weights
echo "[analyze] Computing l3 weights for kappa = ${KAPPA}..."
IFS=',' read -ra KAPPA_VALS <<< "$KAPPA"
for k in "${KAPPA_VALS[@]}"; do
    echo "[analyze]   kappa = $k"
    python3 scripts/add_l3_weight.py "$LO_ROOT" "$RWGT_ROOT" --l3 "$k"
done

# 5. Run plotting scripts
echo "[analyze] Running plotting scripts..."

# 5a. Weight ratio
python3 scripts/plot_weight_ratio.py "$LO_ROOT" "$RWGT_ROOT" --out output/plots/weight_ratio.png

# 5b. C1 vs pT
python3 scripts/plot_C1_vs_pt.py "$LO_ROOT" "$RWGT_ROOT" --out output/plots/C1_vs_pt.png

# 5c. Compare and C1 (multi-variable)
python3 scripts/compare_and_C1.py "$LO_ROOT" "$RWGT_ROOT" --out output/plots/compare_and_C1.png

# 5d. Compare ROOTs (scatter + hist)
python3 scripts/compare_roots.py "$LO_ROOT" "$RWGT_ROOT"

# 5e. Plot kappa3 variants
# Find a generated *_l3corr_*.root file (excluding nloew) to use as the NLO base.
L3CORR_BASE=""
# Prefer the SM (l3=1.0) file if it exists
for f in output/events_l3corr_1p0.root output/events_l3corr_1p00.root; do
    if [[ -f "$f" ]]; then
        L3CORR_BASE="$f"
        break
    fi
done
# Otherwise pick the first non-nloew l3corr file
if [[ -z "$L3CORR_BASE" ]]; then
    for f in output/events_l3corr_*.root; do
        if [[ -f "$f" && "$f" != *"nloew"* ]]; then
            L3CORR_BASE="$f"
            break
        fi
    done
fi

if [[ -n "$L3CORR_BASE" ]]; then
    python3 scripts/plot_kappa3.py "$L3CORR_BASE" "$LO_ROOT" --feature h_pt
else
    echo "[analyze] Warning: no *_l3corr_*.root file found; skipping plot_kappa3.py"
fi

# 6. Move any remaining plots generated in CWD into output/plots/
# compare_roots.py outputs: compare_scatter.png, compare_hpt_hist.png, compare_zpt_hist.png
# compare_roots.py outputs to CWD
for f in compare_scatter.png compare_hpt_hist.png compare_zpt_hist.png; do
    if [[ -f "$f" ]]; then
        mv "$f" output/plots/
    fi
done

# plot_kappa3.py outputs filenames based on the base file name
for f in events_l3corr_*_nlo_only.png events_l3corr_*_vs_lo.png; do
    if [[ -f "$f" ]]; then
        mv "$f" output/plots/
    fi
done

# Also move any compare_and_C1 variant plots that might have ended up in CWD
for f in compare_and_C1_*.png; do
    if [[ -f "$f" ]]; then
        mv "$f" output/plots/
    fi
done

echo "[analyze] Done. Plots saved to output/plots/"
