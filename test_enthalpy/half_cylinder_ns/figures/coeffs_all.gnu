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
set ylabel 'Pressure coefficient Cp'
unset yrange
set title 'quad - Cp' font 'Times,12'
plot 'outputs/quad_three_wave_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 1 title 'three wave', \
     'outputs/quad_three_wave_enthalpy_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 2 title 'three wave enthalpy', \
     'outputs/quad_multi_point_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 3 title 'multi point', \
     'outputs/quad_multi_point_enthalpy_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 4 title 'multi point enthalpy', \
     'fun3D_cp.csv' using 1:2 smooth unique w l lc black lw 2.5 title 'LAURA', \
     'fun3D_cp.csv' using (-$1):2 smooth unique w l lc black lw 2.5 notitle
set ylabel 'Heating rate C_H'
set yrange [0:0.02]
set title 'quad - C_H' font 'Times,12'
plot 'outputs/quad_three_wave_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 1 title 'three wave', \
     'outputs/quad_three_wave_enthalpy_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 2 title 'three wave enthalpy', \
     'outputs/quad_multi_point_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 3 title 'multi point', \
     'outputs/quad_multi_point_enthalpy_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 4 title 'multi point enthalpy', \
     'fun3D_st.csv' using 1:2 smooth unique w l lc black lw 2.5 title 'LAURA', \
     'fun3D_st.csv' using (-$1):2 smooth unique w l lc black lw 2.5 notitle
set ylabel 'Pressure coefficient Cp'
unset yrange
set title 'tri - Cp' font 'Times,12'
plot 'outputs/tri_three_wave_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 1 title 'three wave', \
     'outputs/tri_three_wave_enthalpy_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 2 title 'three wave enthalpy', \
     'outputs/tri_multi_point_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 3 title 'multi point', \
     'outputs/tri_multi_point_enthalpy_o1/coeffs_export.csv' using 7:8 smooth unique w lp lt 4 title 'multi point enthalpy', \
     'fun3D_cp.csv' using 1:2 smooth unique w l lc black lw 2.5 title 'LAURA', \
     'fun3D_cp.csv' using (-$1):2 smooth unique w l lc black lw 2.5 notitle
set ylabel 'Heating rate C_H'
set yrange [0:0.02]
set title 'tri - C_H' font 'Times,12'
plot 'outputs/tri_three_wave_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 1 title 'three wave', \
     'outputs/tri_three_wave_enthalpy_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 2 title 'three wave enthalpy', \
     'outputs/tri_multi_point_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 3 title 'multi point', \
     'outputs/tri_multi_point_enthalpy_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lt 4 title 'multi point enthalpy', \
     'fun3D_st.csv' using 1:2 smooth unique w l lc black lw 2.5 title 'LAURA', \
     'fun3D_st.csv' using (-$1):2 smooth unique w l lc black lw 2.5 notitle
unset multiplot
