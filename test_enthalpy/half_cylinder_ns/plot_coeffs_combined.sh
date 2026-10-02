#!/bin/bash
# The four wall-coefficient plots in a single image: columns are Cp and heating rate C_H,
# rows are the quad and tri meshes. Four separate panels, not overlaid -- with eight curves
# on one set of axes the two meshes of a given scheme share a colour and cannot be told apart.
# Each panel carries the four schemes plus the LAURA reference.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
mkdir -p figures

SCHEMES="three_wave three_wave_enthalpy multi_point multi_point_enthalpy"

{
cat <<'HDR'
set terminal pdf enhanced font 'Times,11' size 9in,6.4in
set output 'figures/coeffs_all.pdf'
set datafile separator ','
# Wong colour-blind-safe palette, one colour per scheme.
set linetype 1 lc rgb '#E69F00' lw 2 pt 7
set linetype 2 lc rgb '#56B4E9' lw 2 pt 9
set linetype 3 lc rgb '#009E73' lw 2 pt 5
set linetype 4 lc rgb '#D55E00' lw 2 pt 13
set grid
set pointsize 0.35
set xrange [-1.6:1.6]
set xlabel 'Angle {/Symbol q}'
set key bottom center font 'Times,9' spacing 1.0 samplen 3
st(x)=2*x/(1e-3*5000**3)
set multiplot layout 2,2
HDR

for mesh in quad tri; do
  for kind in cp ch; do
    if [ "$kind" = cp ]; then
      echo "set ylabel 'Pressure coefficient Cp'"; echo "unset yrange"
      echo "set title '${mesh} - Cp' font 'Times,12'"
      ref=fun3D_cp.csv; using='7:8'
    else
      echo "set ylabel 'Heating rate C_H'"; echo "set yrange [0:0.02]"
      echo "set title '${mesh} - C_H' font 'Times,12'"
      ref=fun3D_st.csv; using='7:(st($9))'
    fi
    lt=0; first=1
    for s in $SCHEMES; do
      lt=$((lt+1))
      f="outputs/${mesh}_${s}_o1/coeffs_export.csv"
      [ -f "$f" ] || continue
      if [ $first = 1 ]; then printf "plot "; first=0; else printf ", \\\\\n     "; fi
      printf "'%s' using %s smooth unique w lp lt %d title '%s'" \
             "$f" "$using" "$lt" "$(echo $s | tr '_' ' ')"
    done
    printf ", \\\\\n     '%s' using 1:2 smooth unique w l lc black lw 2.5 title 'LAURA'" "$ref"
    printf ", \\\\\n     '%s' using (-\$1):2 smooth unique w l lc black lw 2.5 notitle\n" "$ref"
  done
done
echo "unset multiplot"
} > figures/coeffs_all.gnu

gnuplot figures/coeffs_all.gnu && echo "wrote figures/coeffs_all.pdf"
