#!/bin/bash
# setup_env.sh — Detect or install external dependencies for VHtrilinear pipeline
# Safe to re-run (idempotent). Run once before using run_pipeline.sh.
set -euo pipefail

BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$BASEDIR/stage3_fixes.sh"
echo "============================================="
echo " VHtrilinear Environment Setup"
echo " Working directory: $BASEDIR"
echo "============================================="

# ── Helper function ──
check_ok()   { echo "  [✓] $1"; }
check_skip() { echo "  [–] $1 (already present)"; }
check_fail() { echo "  [✗] $1"; exit 1; }

# ── 1. Apptainer ──
echo ""
echo "── Checking Apptainer ──"
if command -v apptainer &>/dev/null; then
    check_ok "apptainer found: $(command -v apptainer)"
else
    check_fail "apptainer not found on PATH. Please install Apptainer first (https://apptainer.org/docs/admin/main/installation.html)"
fi

# ── 2. Container image ──
echo ""
echo "── Checking container image ──"
if [ -f "$BASEDIR/trilinear-boost.sif" ]; then
    check_skip "trilinear-boost.sif"
else
    check_fail "trilinear-boost.sif not found. Build it locally with: apptainer build trilinear-boost.sif docker://lordcristo/trilinear-boost"
fi

# ── 3. LHAPDF 6.2.1 ──
echo ""
echo "── Checking LHAPDF 6.2.1 ──"
LHAPDF_PREFIX="$BASEDIR/lhapdf"
if [ -x "$LHAPDF_PREFIX/bin/lhapdf-config" ]; then
    check_skip "LHAPDF ($LHAPDF_PREFIX/bin/lhapdf-config)"
else
    LHAPDF_VERSION="LHAPDF-6.2.1"
    if [ ! -d "$BASEDIR/$LHAPDF_VERSION" ]; then
        echo "  Downloading $LHAPDF_VERSION..."
        wget -q "https://lhapdf.hepforge.org/downloads/?f=${LHAPDF_VERSION}.tar.gz" -O "$BASEDIR/${LHAPDF_VERSION}.tar.gz"
        tar xf "$BASEDIR/${LHAPDF_VERSION}.tar.gz" -C "$BASEDIR"
        rm -f "$BASEDIR/${LHAPDF_VERSION}.tar.gz"
    fi
    echo "  Building LHAPDF (prefix=$LHAPDF_PREFIX)..."
    mkdir -p "$LHAPDF_PREFIX"
    pushd "$BASEDIR/$LHAPDF_VERSION" > /dev/null
    ./configure --prefix="$LHAPDF_PREFIX" > /dev/null 2>&1
    make -j8 > /dev/null 2>&1
    make install > /dev/null 2>&1
    popd > /dev/null
    check_ok "LHAPDF installed to $LHAPDF_PREFIX"
fi

# Fix LHAPDF URL base for PDF set downloads
LHAPDF_SCRIPT="$LHAPDF_PREFIX/bin/lhapdf"
if [ -f "$LHAPDF_SCRIPT" ]; then
    sed -i "s|^urlbase = .*|urlbase = 'https://lhapdfsets.web.cern.ch/lhapdfsets/current/'|" "$LHAPDF_SCRIPT"
    check_ok "Updated LHAPDF URL base to CERN mirror"
fi

# ── 4. MG5_aMC_v2_5_5 ──
echo ""
echo "── Checking MG5_aMC_v2_5_5 ──"
MG5DIR="$BASEDIR/MG5_aMC_v2_5_5"
if [ -x "$MG5DIR/bin/mg5_aMC" ]; then
    check_skip "MG5_aMC_v2_5_5"
else
    echo "  Downloading MG5_aMC_v2.5.5..."
    wget -q "https://launchpad.net/mg5amcnlo/lts/2.5.x/+download/MG5_aMC_v2.5.5.tar.gz" -O "$BASEDIR/MG5_aMC_v2.5.5.tar.gz"
    tar xf "$BASEDIR/MG5_aMC_v2.5.5.tar.gz" -C "$BASEDIR"
    rm -f "$BASEDIR/MG5_aMC_v2.5.5.tar.gz"
    check_ok "MG5_aMC_v2_5_5 installed"
fi

stage3_patch_mg5_boost_wrappers "$MG5DIR"
check_ok "Patched MG5 PDF wrapper templates for Boost headers"

# ── 5. Configure MG5: lhapdf-config path ──
echo ""
echo "── Configuring MG5 ──"
MG5_CONF="$MG5DIR/input/mg5_configuration.txt"
LHAPDF_CONFIG="$LHAPDF_PREFIX/bin/lhapdf-config"

# Set lhapdf-config path (uncomment and set if commented out, or replace existing)
if grep -q "^# *lhapdf.*=" "$MG5_CONF" 2>/dev/null; then
    # Uncomment and set
    sed -i "s|^# *lhapdf.*=.*|lhapdf = $LHAPDF_CONFIG|" "$MG5_CONF"
    check_ok "Set lhapdf-config in mg5_configuration.txt"
elif grep -q "^lhapdf.*=" "$MG5_CONF" 2>/dev/null; then
    sed -i "s|^lhapdf.*=.*|lhapdf = $LHAPDF_CONFIG|" "$MG5_CONF"
    check_ok "Updated lhapdf-config in mg5_configuration.txt"
else
    echo "lhapdf = $LHAPDF_CONFIG" >> "$MG5_CONF"
    check_ok "Appended lhapdf-config to mg5_configuration.txt"
fi

# Disable COLLIER
if grep -q "^# *collier" "$MG5_CONF" 2>/dev/null; then
    sed -i "s|^# *collier.*|collier = None|" "$MG5_CONF"
    check_ok "Disabled COLLIER"
elif ! grep -q "^collier" "$MG5_CONF" 2>/dev/null; then
    echo "collier = None" >> "$MG5_CONF"
    check_ok "Added collier = None"
else
    check_skip "COLLIER already configured"
fi

# ── 6. Copy UFO model into MG5 models/ ──
echo ""
echo "── Installing UFO model ──"

# Determine the model directory name from the proc card
MODEL_NAME=$(grep "^import model" "$BASEDIR/cards/proc_mc_hz" | awk '{print $3}')
if [ -z "$MODEL_NAME" ]; then
    MODEL_NAME="hhh-model"
fi

# Find the source model directory (could be hhh-model or hhh-model-new)
SRC_MODEL=""
for candidate in "$BASEDIR/trilinear-RW/$MODEL_NAME" "$BASEDIR/trilinear-RW/hhh-model" "$BASEDIR/trilinear-RW/hhh-model-new"; do
    if [ -d "$candidate" ]; then
        SRC_MODEL="$candidate"
        break
    fi
done

if [ -z "$SRC_MODEL" ]; then
    check_fail "Cannot find UFO model directory in trilinear-RW/"
fi

DEST_MODEL="$MG5DIR/models/$MODEL_NAME"
if [ -d "$DEST_MODEL" ]; then
    check_skip "UFO model '$MODEL_NAME' in MG5 models/"
else
    cp -r "$SRC_MODEL" "$DEST_MODEL"
    check_ok "Copied UFO model '$MODEL_NAME' → MG5 models/"
fi

# ── 7. Install diagram filter ──
echo ""
echo "── Installing diagram filter ──"
LOOP_DIR="$MG5DIR/madgraph/loop"
ORIG_FILTER="$LOOP_DIR/loop_diagram_generation.py"
if [ -f "$ORIG_FILTER" ] && [ ! -f "$LOOP_DIR/loop_diagram_generation_original.py" ]; then
    cp "$ORIG_FILTER" "$LOOP_DIR/loop_diagram_generation_original.py"
    check_ok "Backed up original loop_diagram_generation.py"
fi
cp "$BASEDIR/trilinear-RW/vvh-loop_diagram_generation.py" "$ORIG_FILTER"
check_ok "Installed VVH diagram filter"

# ── 8. Generate helper process for PDF libraries ──
echo ""
echo "── Checking PDF libraries ──"
HHH_LIBS="$MG5DIR/HHH-libs"
if [ -f "$HHH_LIBS/libpdf.a" ] && [ -d "$HHH_LIBS/Pdfdata" ]; then
    check_skip "PDF libraries in HHH-libs/"
else
    echo "  Generating helper process (p p > z) to extract PDF libraries..."
    mkdir -p "$HHH_LIBS"
    cat > "$BASEDIR/_tmp_proc_ppz" <<'PROCEOF'
generate p p > z
output ppz
launch ppz -f
quit
PROCEOF
    pushd "$MG5DIR" > /dev/null
    ./bin/mg5_aMC < "$BASEDIR/_tmp_proc_ppz" > /dev/null 2>&1 || true
    popd > /dev/null
    # Extract libraries
    if [ -d "$MG5DIR/ppz/lib" ]; then
        cp "$MG5DIR/ppz/lib/libpdf.a" "$HHH_LIBS/" 2>/dev/null || true
        cp -r "$MG5DIR/ppz/lib/Pdfdata" "$HHH_LIBS/" 2>/dev/null || true
    fi
    # Get libLHAPDF.a
    cp "$LHAPDF_PREFIX/lib/libLHAPDF.a" "$HHH_LIBS/" 2>/dev/null || true
    # Create PDFsets symlink
    ln -sf "$($LHAPDF_CONFIG --datadir)" "$HHH_LIBS/PDFsets" 2>/dev/null || true
    rm -f "$BASEDIR/_tmp_proc_ppz"
    check_ok "PDF libraries prepared in HHH-libs/"
fi

# ── 9. Write lhapdf_env.sh for container usage ──
echo ""
echo "── Writing lhapdf_env.sh ──"
cat > "$BASEDIR/lhapdf_env.sh" <<ENVEOF
#!/bin/bash
# Source this to set up LHAPDF paths, then optionally enter the container.
export PYTHONPATH="$LHAPDF_PREFIX/lib64/python2.7/site-packages/:\${PYTHONPATH:-}"
export LD_LIBRARY_PATH="$LHAPDF_PREFIX/lib/:\${LD_LIBRARY_PATH:-}"
export PATH="$LHAPDF_PREFIX/bin:\$PATH"
ENVEOF
chmod +x "$BASEDIR/lhapdf_env.sh"
check_ok "lhapdf_env.sh written"

# ── Summary ──
echo ""
echo "============================================="
echo " Setup complete!"
echo "============================================="
echo "  Container:  $BASEDIR/trilinear-boost.sif"
echo "  LHAPDF:     $LHAPDF_PREFIX"
echo "  MG5:        $MG5DIR"
echo "  UFO model:  $DEST_MODEL"
echo ""
echo "  Next step:  ./run_pipeline.sh"
echo "============================================="
