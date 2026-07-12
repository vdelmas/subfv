#!/bin/bash

SCHEMES=$(cat ../schemes.txt)
METHODS=$(cat methods.txt)
BASEDIR=$(pwd)

format_title() {
    echo "$1" | tr '_' ' ' | awk '{
        for (i=1; i<=NF; i++) {
            $i = toupper(substr($i,1,1)) substr($i,2)
        }
        print
    }'
}

plot_case() {
    local DIR=$1
    local TITLE=$2

    local FILE="${DIR}/coeffs_export.csv"
    [ -f "${FILE}" ] || return 0

    local GNU="${DIR}/plot_cp.gnu"
    echo "set terminal pdf enhanced font 'Times,10' size 4in,2.25in" > "${GNU}"
    echo "set output '${DIR}/plot_cp.pdf'" >> "${GNU}"
    echo "set datafile separator ','" >> "${GNU}"
    echo "set grid" >> "${GNU}"
    echo "set pointsize 0.5" >> "${GNU}"
    echo "set key bottom center" >> "${GNU}"
    echo "set xrange [-1.6:1.6]" >> "${GNU}"
    echo "set xlabel 'Angle {/Symbol q}'" >> "${GNU}"
    echo "set ylabel 'Pressure coefficient Cp'" >> "${GNU}"
    echo -n "plot '${FILE}' using 7:8 smooth unique w lp lw 2 title '${TITLE}'" >> "${GNU}"
    echo -n ", '${BASEDIR}/fun3D_cp.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA'" >> "${GNU}"
    echo    ", '${BASEDIR}/fun3D_cp.csv' using (-\$1):2 smooth unique w l lc black lw 2 notitle" >> "${GNU}"
    gnuplot "${GNU}"

    GNU="${DIR}/plot_st.gnu"
    echo "set terminal pdf enhanced font 'Times,10' size 4in,2.25in" > "${GNU}"
    echo "set output '${DIR}/plot_st.pdf'" >> "${GNU}"
    echo "set datafile separator ','" >> "${GNU}"
    echo "set grid" >> "${GNU}"
    echo "set pointsize 0.5" >> "${GNU}"
    echo "set key bottom center" >> "${GNU}"
    echo "set xrange [-1.6:1.6]" >> "${GNU}"
    echo "set yrange [0:0.02]" >> "${GNU}"
    echo "st(x)=2*x/(1e-3*5000**3)" >> "${GNU}"
    echo "set xlabel 'Angle {/Symbol q}'" >> "${GNU}"
    echo "set ylabel 'Heating rate C_H (q/(1/2 {/Symbol r}V^3))'" >> "${GNU}"
    echo -n "plot '${FILE}' using 7:(st(\$9)) smooth unique w lp lw 2 title '${TITLE}'" >> "${GNU}"
    echo -n ", '${BASEDIR}/fun3D_st.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA'" >> "${GNU}"
    echo    ", '${BASEDIR}/fun3D_st.csv' using (-\$1):2 smooth unique w l lc black lw 2 notitle" >> "${GNU}"
    gnuplot "${GNU}"

    echo "  -> ${DIR}/plot_cp.pdf + plot_st.pdf"
}

echo "=== O1 ==="
for SCHEME in ${SCHEMES}; do
    TITLE=$(format_title "${SCHEME}")
    plot_case "outputs/${SCHEME}" "${TITLE}"
done

echo "=== O2 ==="
for METHOD in ${METHODS}; do
    for SCHEME in ${SCHEMES}; do
        TITLE=$(format_title "${SCHEME} O2 ${METHOD}")
        plot_case "outputs/${SCHEME}_o2_${METHOD}" "${TITLE}"
    done
done
