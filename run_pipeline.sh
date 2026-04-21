#!/bin/bash
# run_pipeline.sh — Generate LO ZH events and reweight for λ₃ NLO EW corrections.
# Prerequisites: run setup_env.sh first.
# Usage: ./run_pipeline.sh [--nevents 500000] [--ecm 13600] [--skip-generation]
set -euo pipefail

BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REAL_BASEDIR=$(readlink -f "$BASEDIR")
source "$REAL_BASEDIR/stage3_fixes.sh"

CONTAINER_IMAGE="${TRILINEAR_SIF:-$REAL_BASEDIR/trilinear-boost.sif}"
if [ ! -f "$CONTAINER_IMAGE" ]; then
    echo "[✗] Missing container image: $CONTAINER_IMAGE"
    echo "    Build the Boost-enabled SIF with:"
    echo "    apptainer build trilinear-boost.sif docker://lordcristo/trilinear-boost"
    exit 1
fi

# ── Enter Apptainer container if not already inside ──
if [ -z "${APPTAINER_CONTAINER:-}" ] && [ -z "${SINGULARITY_CONTAINER:-}" ]; then
    echo "── Entering Apptainer container ──"
    # Launch from the real parent directory to avoid AFS symlink issues
    (cd "$REAL_BASEDIR/.." && \
     apptainer exec --env PATH="$PATH" --env LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" --env PYTHONPATH="${PYTHONPATH:-}" \
     "$CONTAINER_IMAGE" bash "$REAL_BASEDIR/run_pipeline.sh" "$@")
    exit $?
fi

cd "$BASEDIR"

source "$BASEDIR/lhapdf_env.sh"
stage3_verify_container_boost_headers

# ── Parse arguments ──
NEVENTS=500000
ECM=13600       # in GeV (13.6 TeV = 13600)
SKIP_GEN=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --nevents)    NEVENTS="$2"; shift 2 ;;
        --ecm)        ECM="$2"; shift 2 ;;
        --skip-generation) SKIP_GEN=true; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

# Beam energy per proton = ECM / 2
EBEAM=$(echo "$ECM / 2" | bc)

echo "============================================="
echo " VHtrilinear Run Pipeline"
echo " Events: $NEVENTS    √s: $ECM GeV"
echo "============================================="

MG5DIR="$BASEDIR/MG5_aMC_v2_5_5"
RWDIR="$BASEDIR/trilinear-RW"
LHAPDF_CONFIG="$BASEDIR/lhapdf/bin/lhapdf-config"

# Determine model name from proc card
MODEL_NAME=$(grep "^import model" "$BASEDIR/cards/proc_mc_hz" | awk '{print $3}')

# ── Step 1: Generate hz_MC (LO QCD process) ──
echo ""
echo "── Step 1: Generating hz_MC ──"
stage3_patch_mg5_boost_wrappers "$MG5DIR"
if [ -d "$MG5DIR/hz_MC" ] && [ "$SKIP_GEN" = true ]; then
    echo "  [–] hz_MC already exists (skipping)"
else
    rm -rf "$MG5DIR/hz_MC"
    pushd "$MG5DIR" > /dev/null
    ./bin/mg5_aMC < "$BASEDIR/cards/proc_mc_hz"
    popd > /dev/null
    echo "  [✓] hz_MC generated"
fi
stage3_patch_mg5_boost_wrappers "$MG5DIR"

# Keep the generated MG5 config in sync with the local LHAPDF install.
HZ_MC_AMCNLO_CONF="$MG5DIR/hz_MC/Cards/amcatnlo_configuration.txt"
if [ -f "$HZ_MC_AMCNLO_CONF" ]; then
    sed -i "s|^lhapdf = .*|lhapdf = $LHAPDF_CONFIG #|" "$HZ_MC_AMCNLO_CONF"
    echo "  [✓] amcatnlo_configuration.txt patched with LHAPDF path"
fi

# ── Step 2: Configure run_card.dat ──
echo ""
echo "── Step 2: Configuring run_card.dat ──"
RUNCARD="$MG5DIR/hz_MC/Cards/run_card.dat"

# Copy our template run_card over the default
cp "$BASEDIR/cards/run_card.dat" "$RUNCARD"

# Apply dynamic settings
sed -i "s|.*=.*nevents.*|  $NEVENTS = nevents ! Number of unweighted events requested|" "$RUNCARD"
sed -i "s|.*=.*ebeam1.*|  $EBEAM = ebeam1  ! beam 1 total energy in GeV|" "$RUNCARD"
sed -i "s|.*=.*ebeam2.*|  $EBEAM = ebeam2  ! beam 2 total energy in GeV|" "$RUNCARD"

# Ensure store_rwgt_info is True
sed -i 's|False.*=.*store_rwgt_inf|True = store_rwgt_inf|' "$RUNCARD"

echo "  [✓] run_card.dat configured ($NEVENTS events, √s=$ECM GeV)"

# Copy param_card
cp "$BASEDIR/cards/param_card.dat" "$MG5DIR/hz_MC/Cards/param_card.dat"
echo "  [✓] param_card.dat copied"

# ── Step 3: Run gevirt.sh ──
echo ""
echo "── Step 3: Extracting virtual EW subprocesses ──"
cp "$RWDIR/gevirt.sh" "$MG5DIR/"
pushd "$MG5DIR" > /dev/null
./gevirt.sh hz_MC
popd > /dev/null

if [ ! -f "$MG5DIR/check_olp.inc" ] || [ ! -f "$MG5DIR/proc_ml" ]; then
    echo "  [✗] gevirt.sh failed — check_olp.inc or proc_ml not found"
    exit 1
fi
echo "  [✓] check_olp.inc and proc_ml generated"

# ── Step 4: Generate hz_ME (EW virtual matrix elements) ──
echo ""
echo "── Step 4: Generating hz_ME ──"
rm -rf "$MG5DIR/hz_ME"

# Build proc card for ME generation
cat > "$MG5DIR/_proc_hz_me" <<EOF
import model $MODEL_NAME
$(cat "$MG5DIR/proc_ml")
output hz_ME
quit
EOF

pushd "$MG5DIR" > /dev/null
./bin/mg5_aMC < _proc_hz_me
rm -f _proc_hz_me
popd > /dev/null
echo "  [✓] hz_ME generated"
stage3_patch_mg5_boost_wrappers "$MG5DIR"

# ── Step 5: Build check_OLP executable ──
echo ""
echo "── Step 5: Building reweighting executable ──"
SUBPROC="$MG5DIR/hz_ME/SubProcesses"

# Copy reweighting code
cp "$RWDIR/makefile"    "$SUBPROC/"
cp "$RWDIR/check_OLP.f" "$SUBPROC/"
stage3_patch_check_olp_pdflabel "$SUBPROC/check_OLP.f"

# Copy generated include
cp "$MG5DIR/check_olp.inc" "$SUBPROC/"

# Auto-detect subprocess directory name (e.g. P0_uux_hz)
ME_SUBDIR=$(find "$SUBPROC" -maxdepth 1 -type d -name 'P0_*' | head -1)
if [ -z "$ME_SUBDIR" ]; then
    echo "  [✗] No P0_* subprocess directory found in hz_ME/SubProcesses/"
    exit 1
fi
ME_SUBDIR_NAME=$(basename "$ME_SUBDIR")
echo "  Using ME subprocess: $ME_SUBDIR_NAME"

cp "$ME_SUBDIR/pmass.inc"        "$SUBPROC/"
cp "$ME_SUBDIR/nsqso_born.inc"   "$SUBPROC/"
cp "$ME_SUBDIR/nsquaredSO.inc"   "$SUBPROC/"

# Auto-detect MC subprocess directory
MC_SUBDIR=$(find "$MG5DIR/hz_MC/SubProcesses" -maxdepth 1 -type d -name 'P0_*' | head -1)
if [ -z "$MC_SUBDIR" ]; then
    echo "  [✗] No P0_* subprocess directory found in hz_MC/SubProcesses/"
    exit 1
fi

cp "$MG5DIR/hz_MC/SubProcesses/c_weight.inc" "$SUBPROC/"
cp "$MC_SUBDIR/nexternal.inc"               "$SUBPROC/"

# PDF libraries
LIBS="$MG5DIR/hz_ME/lib"
cp "$MG5DIR/HHH-libs/libpdf.a"    "$LIBS/" 2>/dev/null || true
cp -r "$MG5DIR/HHH-libs/Pdfdata"  "$LIBS/" 2>/dev/null || true
cp "$MG5DIR/HHH-libs/libLHAPDF.a" "$LIBS/" 2>/dev/null || true
if [ -L "$MG5DIR/HHH-libs/PDFsets" ] || [ -d "$MG5DIR/HHH-libs/PDFsets" ]; then
    cp -r "$MG5DIR/HHH-libs/PDFsets" "$LIBS/" 2>/dev/null || true
fi

# Build
pushd "$SUBPROC" > /dev/null
make OLP_static
make check_OLP
popd > /dev/null

if [ ! -x "$SUBPROC/check_OLP" ]; then
    echo "  [✗] check_OLP build failed"
    exit 1
fi
echo "  [✓] check_OLP executable built"

# ── Step 6: Generate LO events ──
echo ""
echo "── Step 6: Generating LO events ──"
cat > "$MG5DIR/_launch_hz" <<'LAUNCHEOF'
launch hz_MC
noshowerLO
done
quit
LAUNCHEOF

pushd "$MG5DIR" > /dev/null
./bin/mg5_aMC < _launch_hz
rm -f _launch_hz
popd > /dev/null

# Find the generated LHE file
LHE_GZ=$(find "$MG5DIR/hz_MC/Events/" -name "events.lhe.gz" -o -name "unweighted_events.lhe.gz" | head -1)
if [ -z "$LHE_GZ" ]; then
    echo "  [✗] No events.lhe.gz found in hz_MC/Events/"
    exit 1
fi
echo "  [✓] LO events generated: $LHE_GZ"

# ── Step 7: Reweight ──
echo ""
echo "── Step 7: Reweighting events ──"
gunzip -c "$LHE_GZ" > "$SUBPROC/events.lhe"
echo "  Unzipped $(wc -l < "$SUBPROC/events.lhe") lines to events.lhe"

pushd "$SUBPROC" > /dev/null
./check_OLP
popd > /dev/null

if [ ! -f "$SUBPROC/events_rwgt.lhe" ]; then
    echo "  [✗] Reweighting failed — events_rwgt.lhe not found"
    exit 1
fi
echo "  [✓] Reweighting complete"

# ── Step 8: Copy to output/ ──
echo ""
echo "── Copying results to output/ ──"
mkdir -p "$BASEDIR/output"
cp "$SUBPROC/events.lhe"      "$BASEDIR/output/events.lhe"
cp "$SUBPROC/events_rwgt.lhe" "$BASEDIR/output/events_rwgt.lhe"
echo "  [✓] events.lhe → output/"
echo "  [✓] events_rwgt.lhe → output/"

echo ""
echo "============================================="
echo " Pipeline complete!"
echo " Output: $BASEDIR/output/"
echo " Next:   ./analyze.sh"
echo "============================================="
