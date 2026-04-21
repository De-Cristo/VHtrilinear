#!/bin/bash
# analyze.sh — Stage 4 wrapper: runs the unified Python analysis pipeline.
# Delegates all work to scripts/analyze.py.
#
# Usage: ./analyze.sh [--kappa 0,1,2,5,10,-2,-5,-10] [--ebeam 6800] [--lo-lhe ...] [--rw-lhe ...]
set -euo pipefail

BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Pass all arguments straight through to analyze.py
python3 "$BASEDIR/scripts/analyze.py" \
    --lo-lhe "$BASEDIR/output/events.lhe" \
    --rw-lhe "$BASEDIR/output/events_rwgt.lhe" \
    --outdir "$BASEDIR/output" \
    "$@"
