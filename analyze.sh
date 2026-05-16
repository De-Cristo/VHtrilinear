#!/bin/bash
# analyze.sh — Stage 4 wrapper: runs the unified Python analysis pipeline.
# Delegates all work to scripts/analyze.py.
#
# Usage: ./analyze.sh [--process zh|wh] [--kappa 0,1,2,5,10,-2,-5,-10] [--ebeam 6800]
set -euo pipefail

BASEDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROCESS=zh
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --process) PROCESS="$2"; shift 2 ;;
        *) ARGS+=("$1"); shift ;;
    esac
done

python3 "$BASEDIR/scripts/analyze.py" \
    --process "$PROCESS" \
    "${ARGS[@]}"
