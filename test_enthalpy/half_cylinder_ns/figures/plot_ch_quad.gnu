set terminal pdf enhanced font 'Times,12' size 4in,3in
set linetype 1 lc rgb '#E69F00' lw 2 pt 7
set linetype 2 lc rgb '#56B4E9' lw 2 pt 9
set linetype 3 lc rgb '#009E73' lw 2 pt 5
set linetype 4 lc rgb '#D55E00' lw 2 pt 13
set output 'figures/plot_ch_quad.pdf'
set datafile separator ','
set grid
set pointsize 0.4
set key bottom center
set xrange [-1.6:1.6]
set xlabel 'Angle {/Symbol q}'
set ylabel 'Heating rate C_H'
set yrange [0:0.02]
st(x)=2*x/(1e-3*5000**3)
plot 'outputs/quad_three_wave_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'three wave', \
     'outputs/quad_three_wave_enthalpy_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'three wave enthalpy', \
     'outputs/quad_multi_point_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'multi point', \
     'outputs/quad_multi_point_enthalpy_o1/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'multi point enthalpy', \
     'fun3D_st.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA', \
     'fun3D_st.csv' using (-$1):2 smooth unique w l lc black lw 2 notitle
