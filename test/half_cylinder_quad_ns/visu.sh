#!/bin/bash
set -e

SCHEME=$1
NPROC=$2
SECOND=$3
METHOD=$4

BASEDIR=$(pwd)

if [[ "$SECOND" == ".true." ]]; then
  OUTDIR=outputs/${SCHEME}_o2_${METHOD}
else
  OUTDIR=outputs/${SCHEME}
fi

mkdir -p $OUTDIR
cd $OUTDIR

pvbatch ../../visualize_fields.py

if [ -f coeffs_export.csv ]; then
    TITLE=$(echo "${SCHEME}" | tr '_' ' ')

    cat > plot_cp.gnu <<EOF
set terminal pdf enhanced font 'Times,10' size 4in,2.25in
set output 'plot_cp.pdf'
set datafile separator ','
set grid
set pointsize 0.5
set key bottom center
set xrange [-1.6:1.6]
set xlabel 'Angle {/Symbol q}'
set ylabel 'Pressure coefficient Cp'
plot 'coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title '${TITLE}', \\
     '${BASEDIR}/fun3D_cp.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA', \\
     '${BASEDIR}/fun3D_cp.csv' using (-\$1):2 smooth unique w l lc black lw 2 notitle
EOF
    gnuplot plot_cp.gnu

    cat > plot_st.gnu <<EOF
set terminal pdf enhanced font 'Times,10' size 4in,2.25in
set output 'plot_st.pdf'
set datafile separator ','
set grid
set pointsize 0.5
set key bottom center
set xrange [-1.6:1.6]
set yrange [0:0.02]
st(x)=2*x/(1e-3*5000**3)
set xlabel 'Angle {/Symbol q}'
set ylabel 'Heating rate C_H (q/(1/2 {/Symbol r}V^3))'
plot 'coeffs_export.csv' using 7:(st(\$9)) smooth unique w lp lw 2 title '${TITLE}', \\
     '${BASEDIR}/fun3D_st.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA', \\
     '${BASEDIR}/fun3D_st.csv' using (-\$1):2 smooth unique w l lc black lw 2 notitle
EOF
    gnuplot plot_st.gnu
fi

date > visu.timestamp
