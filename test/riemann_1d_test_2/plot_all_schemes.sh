#!/bin/bash
# Generates one PDF per scheme (2x2 multiplot: density, velocity, pressure, ie)
# Usage: bash plot_all_schemes.sh [SCHEME]
#   With SCHEME: generate only that scheme's PDF.
#   Without SCHEME: generate PDFs for all schemes in schemes.txt.
# Output: outputs/<SCHEME>/plot_toro.pdf
set -euo pipefail

if [ $# -ge 1 ]; then
    SCHEMES="$1"
else
    SCHEMES=$(cat ../schemes.txt)
fi

ie_formula="ie(rho,p) = p/(rho*0.4)"

for SCHEME in ${SCHEMES}; do
    FILE=$(find "outputs/${SCHEME}" -name "output_-1.dat" 2>/dev/null | head -1)
    [ -f "${FILE}" ] || continue

    OUTDIR=$(dirname "${FILE}")
    GNU="${OUTDIR}/plot_toro.gnu"
    TITLE=$(echo "${SCHEME}" | tr '_' ' ')

    cat > "${GNU}" <<GNU
set terminal pdf enhanced font 'Times,12' size 8in,6in
set output '${OUTDIR}/plot_toro.pdf'
set multiplot layout 2,2 title '${TITLE}'
${ie_formula}
set grid
set pointsize 0.3
set xlabel 'x'

set ylabel 'Density'
set key top right
plot '${FILE}' using 1:9 smooth unique w lp lw 2 notitle

set ylabel 'Velocity-X'
set key bottom left
plot '${FILE}' using 1:10 smooth unique w lp lw 2 notitle

set ylabel 'Pressure'
set key top right
plot '${FILE}' using 1:13 smooth unique w lp lw 2 notitle

set ylabel 'Internal energy'
set key top left
plot '${FILE}' using 1:(ie(\$9,\$13)) smooth unique w lp lw 2 notitle

unset multiplot
GNU

    gnuplot "${GNU}"
    echo "Generated: ${OUTDIR}/plot_toro.pdf"
done
