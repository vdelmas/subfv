#!/bin/bash
# Cp and heating-rate C_H along the cylinder wall, every scheme against the LAURA reference
# (fun3D_cp.csv / fun3D_st.csv), for one mesh kind.
#
# usage: ./plot_coeffs.sh <quad|tri>
set -e
MESH=$1
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

SCHEMES="three_wave three_wave_enthalpy multi_point multi_point_enthalpy"
TITLES="three_wave three_wave_enthalpy multi_point multi_point_enthalpy"

gen() {   # $1 = cp|ch
  local kind=$1 out gnu ylab ref using
  gnu="figures/plot_${kind}_${MESH}.gnu"
  out="figures/plot_${kind}_${MESH}.pdf"
  if [ "$kind" = "cp" ]; then
    ylab="Pressure coefficient Cp"; ref="fun3D_cp.csv"; using="7:8"
  else
    ylab="Heating rate C_H"; ref="fun3D_st.csv"; using="7:(st(\$9))"
  fi
  {
    echo "set terminal pdf enhanced font 'Times,12' size 4in,3in"
    # Wong colour-blind-safe palette, the project's standard for report plots.
    echo "set linetype 1 lc rgb '#E69F00' lw 2 pt 7"
    echo "set linetype 2 lc rgb '#56B4E9' lw 2 pt 9"
    echo "set linetype 3 lc rgb '#009E73' lw 2 pt 5"
    echo "set linetype 4 lc rgb '#D55E00' lw 2 pt 13"
    echo "set output '$out'"
    echo "set datafile separator ','"
    echo "set grid"
    echo "set pointsize 0.4"
    echo "set key bottom center"
    echo "set xrange [-1.6:1.6]"
    echo "set xlabel 'Angle {/Symbol q}'"
    echo "set ylabel '$ylab'"
    [ "$kind" = "ch" ] && echo "set yrange [0:0.02]"
    echo "st(x)=2*x/(1e-3*5000**3)"
    first=1
    for s in $SCHEMES; do
      f="outputs/${MESH}_${s}_o1/coeffs_export.csv"
      [ -f "$f" ] || continue
      if [ $first = 1 ]; then printf "plot "; first=0; else printf ", \\\\\n     "; fi
      printf "'%s' using %s smooth unique w lp lw 2 title '%s'" "$f" "$using" "$(echo $s | tr '_' ' ')"
    done
    [ $first = 1 ] && { echo "# no data"; return; }
    printf ", \\\\\n     '%s' using 1:2 smooth unique w l lc black lw 2 title 'LAURA'" "$ref"
    printf ", \\\\\n     '%s' using (-\$1):2 smooth unique w l lc black lw 2 notitle\n" "$ref"
  } > "$gnu"
  grep -q '^plot' "$gnu" && gnuplot "$gnu" && echo "wrote $out"
}

mkdir -p figures
gen cp
gen ch
