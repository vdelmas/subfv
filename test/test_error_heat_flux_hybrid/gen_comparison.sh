#!/bin/bash
# Generate comparison.pdf from existing simulation outputs.
# Assumes visu.sh has already been run for each scheme (outputs/SCHEME/*/error.txt exist).
# Usage: bash gen_comparison.sh [SCHEME...]
#   No args: process all schemes listed in schemes.txt
#   With args: process only the given scheme(s)
set -euo pipefail

ROOT=$(pwd)

if [ $# -ge 1 ]; then
    SCHEMES="$*"
else
    SCHEMES=$(grep -v '^#\|^$' schemes.txt)
fi

for SCHEME in $SCHEMES; do
    echo "=== $SCHEME ==="

    bash collect_errors.sh    "$SCHEME"
    bash collect_residuals.sh "$SCHEME"

    gnuplot -e "SCHEME='${SCHEME}'" plot_errors.gnu
    echo "  -> residuals_${SCHEME}.pdf  errors_cp_${SCHEME}.pdf  errors_ch_${SCHEME}.pdf"
done

# Generate corner mesh images if not already done
NEED_CORNERS=0
for label in top_left top_right bottom_left bottom_right; do
    [ -f "corners_visu/${label}_mesh.png" ] || NEED_CORNERS=1
done
if [ $NEED_CORNERS -eq 1 ]; then
    echo "=== Generating corner mesh images ==="
    pvbatch visu_corners.py
else
    echo "=== Corner mesh images already present, skipping ==="
fi

echo "=== Compiling comparison.tex ==="
latexmk -pdf -interaction=nonstopmode comparison.tex
echo "=== Done: comparison.pdf ==="
