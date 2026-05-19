#!/bin/bash
# run_pipeline.sh — Generate LO VH events and reweight for λ₃ NLO EW corrections.
# Prerequisites: run setup_env.sh first.
#
# Usage:
#   ./run_pipeline.sh --process zh [--nevents 500000] [--ecm 13600]
#   ./run_pipeline.sh --process wh [--nevents 500000] [--ecm 13600]
#
# zh:
#   - runs one subchannel
#   - publishes output/zh/
#
# wh:
#   - runs wh_plus and wh_minus internally
#   - publishes only output/wh/
#
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
PROCESS="zh"
NEVENTS=500000
ECM=13600
SKIP_GEN=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        --process) PROCESS="$2"; shift 2 ;;
        --nevents) NEVENTS="$2"; shift 2 ;;
        --ecm) ECM="$2"; shift 2 ;;
        --skip-generation) SKIP_GEN=true; shift ;;
        *) echo "Unknown arg: $1"; exit 1 ;;
    esac
done

process_python() {
    python3 - "$@"
}

get_process_value() {
    local expr="$1"
    process_python "$PROCESS" "$expr" <<'PY'
from pathlib import Path
import sys
from scripts.vh_processes import get_public_process

process_key = sys.argv[1]
expr = sys.argv[2]
spec = get_public_process(process_key)
print(eval(expr))
PY
}

# Beam energy per proton = ECM / 2
EBEAM=$(echo "$ECM / 2" | bc)

echo "============================================="
echo " VHtrilinear Run Pipeline"
echo " Process: $PROCESS    Events: $NEVENTS    √s: $ECM GeV"
echo "============================================="

MG5DIR="$BASEDIR/MG5_aMC_v2_5_5"
RWDIR="$BASEDIR/trilinear-RW"
LHAPDF_CONFIG="$BASEDIR/lhapdf/bin/lhapdf-config"

MODEL_NAME="hhh-model-new"

mapfile -t SUBCHANNEL_KEYS < <(python3 - "$PROCESS" <<'PY'
import sys
from scripts.vh_processes import get_public_process

for sub in get_public_process(sys.argv[1]).subchannels:
    print(sub.key)
PY
)

for SUBCHANNEL in "${SUBCHANNEL_KEYS[@]}"; do
    SUB_MC_DIR=$(python3 - "$PROCESS" "$SUBCHANNEL" <<'PY'
import sys
from scripts.vh_processes import get_public_process

process_key, sub_key = sys.argv[1], sys.argv[2]
for sub in get_public_process(process_key).subchannels:
    if sub.key == sub_key:
        print(sub.mc_dir)
        break
PY
)

    SUB_ME_DIR=$(python3 - "$PROCESS" "$SUBCHANNEL" <<'PY'
import sys
from scripts.vh_processes import get_public_process

process_key, sub_key = sys.argv[1], sys.argv[2]
for sub in get_public_process(process_key).subchannels:
    if sub.key == sub_key:
        print(sub.me_dir)
        break
PY
)

    PROC_CARD=$(python3 - "$PROCESS" "$SUBCHANNEL" <<'PY'
import sys
from scripts.vh_processes import get_public_process

process_key, sub_key = sys.argv[1], sys.argv[2]
repo = get_public_process(process_key)
for sub in repo.subchannels:
    if sub.key == sub_key:
        print(sub.proc_card)
        break
PY
)

    SUB_RUN_CARD=$(python3 - "$PROCESS" <<'PY'
import sys
from scripts.vh_processes import get_public_process

print(get_public_process(sys.argv[1]).run_card)
PY
)

    SUB_PARAM_CARD=$(python3 - "$PROCESS" <<'PY'
import sys
from scripts.vh_processes import get_public_process

print(get_public_process(sys.argv[1]).param_card)
PY
)

    rm -rf "$MG5DIR/$SUB_MC_DIR" "$MG5DIR/$SUB_ME_DIR"

    stage3_patch_mg5_boost_wrappers "$MG5DIR"
    pushd "$MG5DIR" > /dev/null
    ./bin/mg5_aMC < "$BASEDIR/$PROC_CARD"
    popd > /dev/null
    stage3_patch_mg5_boost_wrappers "$MG5DIR"

    SUB_AMCNLO_CONF="$MG5DIR/$SUB_MC_DIR/Cards/amcatnlo_configuration.txt"
    if [ -f "$SUB_AMCNLO_CONF" ]; then
        sed -i "s|^lhapdf = .*|lhapdf = $LHAPDF_CONFIG #|" "$SUB_AMCNLO_CONF"
    fi

    RUNCARD="$MG5DIR/$SUB_MC_DIR/Cards/run_card.dat"
    cp "$BASEDIR/$SUB_RUN_CARD" "$MG5DIR/$SUB_MC_DIR/Cards/run_card.dat"
    cp "$BASEDIR/$SUB_PARAM_CARD" "$MG5DIR/$SUB_MC_DIR/Cards/param_card.dat"
    sed -i "s|.*=.*nevents.*|  $NEVENTS = nevents ! Number of unweighted events requested|" "$RUNCARD"
    sed -i "s|.*=.*ebeam1.*|  $EBEAM = ebeam1  ! beam 1 total energy in GeV|" "$RUNCARD"
    sed -i "s|.*=.*ebeam2.*|  $EBEAM = ebeam2  ! beam 2 total energy in GeV|" "$RUNCARD"
    sed -i 's|False.*=.*store_rwgt_inf|True = store_rwgt_inf|' "$RUNCARD"

    cp "$RWDIR/gevirt.sh" "$MG5DIR/"
    pushd "$MG5DIR" > /dev/null
    ./gevirt.sh "$SUB_MC_DIR"
    popd > /dev/null

    cat > "$MG5DIR/_proc_${SUBCHANNEL}_me" <<EOF
import model $MODEL_NAME
$(cat "$MG5DIR/proc_ml")
output $SUB_ME_DIR
quit
EOF

    pushd "$MG5DIR" > /dev/null
    ./bin/mg5_aMC < "_proc_${SUBCHANNEL}_me"
    rm -f "_proc_${SUBCHANNEL}_me"
    popd > /dev/null
    stage3_patch_mg5_boost_wrappers "$MG5DIR"

    SUBPROC="$MG5DIR/$SUB_ME_DIR/SubProcesses"
    cp "$RWDIR/makefile" "$SUBPROC/"
    cp "$RWDIR/check_OLP.f" "$SUBPROC/"
    stage3_patch_check_olp_pdflabel "$SUBPROC/check_OLP.f"
    cp "$MG5DIR/check_olp.inc" "$SUBPROC/"

    ME_SUBDIR=$(find "$SUBPROC" -maxdepth 1 -type d -name 'P0_*' | head -1)
    MC_SUBDIR=$(find "$MG5DIR/$SUB_MC_DIR/SubProcesses" -maxdepth 1 -type d -name 'P0_*' | head -1)

    cp "$ME_SUBDIR/pmass.inc" "$SUBPROC/"
    cp "$ME_SUBDIR/nsqso_born.inc" "$SUBPROC/"
    cp "$ME_SUBDIR/nsquaredSO.inc" "$SUBPROC/"
    cp "$MG5DIR/$SUB_MC_DIR/SubProcesses/c_weight.inc" "$SUBPROC/"
    cp "$MC_SUBDIR/nexternal.inc" "$SUBPROC/"

    LIBS="$MG5DIR/$SUB_ME_DIR/lib"
    cp "$MG5DIR/HHH-libs/libpdf.a" "$LIBS/" 2>/dev/null || true
    cp -r "$MG5DIR/HHH-libs/Pdfdata" "$LIBS/" 2>/dev/null || true
    cp "$MG5DIR/HHH-libs/libLHAPDF.a" "$LIBS/" 2>/dev/null || true
    if [ -L "$MG5DIR/HHH-libs/PDFsets" ] || [ -d "$MG5DIR/HHH-libs/PDFsets" ]; then
        cp -r "$MG5DIR/HHH-libs/PDFsets" "$LIBS/" 2>/dev/null || true
    fi

    pushd "$SUBPROC" > /dev/null
    make OLP_static
    make check_OLP
    popd > /dev/null

    cat > "$MG5DIR/_launch_${SUBCHANNEL}" <<EOF
launch $SUB_MC_DIR
noshowerLO
done
quit
EOF

    pushd "$MG5DIR" > /dev/null
    ./bin/mg5_aMC < "_launch_${SUBCHANNEL}"
    rm -f "_launch_${SUBCHANNEL}"
    popd > /dev/null

    LHE_GZ=$(find "$MG5DIR/$SUB_MC_DIR/Events/" -name "events.lhe.gz" -o -name "unweighted_events.lhe.gz" | head -1)
    gunzip -c "$LHE_GZ" > "$SUBPROC/events.lhe"

    pushd "$SUBPROC" > /dev/null
    ./check_OLP
    popd > /dev/null

    SUB_OUTDIR="$BASEDIR/output/_${PROCESS}_internal/${SUBCHANNEL}"
    mkdir -p "$SUB_OUTDIR"
    cp "$SUBPROC/events.lhe" "$SUB_OUTDIR/events.lhe"
    cp "$SUBPROC/events_rwgt.lhe" "$SUB_OUTDIR/events_rwgt.lhe"
done

if [[ "$PROCESS" == "zh" ]]; then
    mkdir -p "$BASEDIR/output/zh"
    cp "$BASEDIR/output/_zh_internal/zh/events.lhe" "$BASEDIR/output/zh/events.lhe"
    cp "$BASEDIR/output/_zh_internal/zh/events_rwgt.lhe" "$BASEDIR/output/zh/events_rwgt.lhe"
fi

echo ""
echo "============================================="
echo " Pipeline complete!"
echo " Output: $BASEDIR/output/"
echo " Next:   ./analyze.sh --process $PROCESS"
echo "============================================="
