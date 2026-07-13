set terminal pdf enhanced font 'Times,12' size 4in,3in
set linetype 1 lc rgb "#E69F00" lw 2 pt 7
set linetype 2 lc rgb "#56B4E9" lw 2 pt 9
set linetype 3 lc rgb "#009E73" lw 2 pt 5
set linetype 4 lc rgb "#F0E442" lw 2 pt 11
set linetype 5 lc rgb "#0072B2" lw 2 pt 13
set linetype 6 lc rgb "#D55E00" lw 2 pt 6
set linetype 7 lc rgb "#CC79A7" lw 2 pt 8
set linetype 8 lc rgb "#000000" lw 2 pt 4
set output 'plot_cp.pdf'
set datafile separator ','
set grid
set xlabel 'Angle {/Symbol q}'
set ylabel 'Pressure coefficient Cp'
set xrange [-1.6:1.6]
set key bottom center
set pointsize 0.5
plot 'outputs/two_wave/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'Two Wave', 'outputs/three_wave/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'Three Wave', 'outputs/modified_three_wave/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'Modified Three Wave', 'outputs/multi_point/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'Multi Point', 'outputs/multi_point_iso/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'Multi Point Iso', 'outputs/ZB_AR1D_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB AR1D LS', 'outputs/ZB_ARMD_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB ARMD LS', 'outputs/ZB_ARMD_LSM/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB ARMD LSM', 'outputs/ZB_ARMDM_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB ARMDM LS', 'outputs/ZB_ARMDMAT_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB ARMDMAT LS', 'outputs/ZB_ARMDMMAT_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB ARMDMMAT LS', 'outputs/ZB_AM_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB AM LS', 'outputs/ZB_AM_LSM/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB AM LSM', 'outputs/ZB_AMISO_LS/coeffs_export.csv' using 7:8 smooth unique w lp lw 2 title 'ZB AMISO LS','fun3D_cp.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA','fun3D_cp.csv' using (-$1):2 smooth unique w l lc black lw 2 notitle
set output 'plot_q.pdf'
set datafile separator ','
set grid
set xlabel 'Angle {/Symbol q}'
set ylabel 'Heat flux (W/m^2)'
set xrange [-1.6:1.6]
set yrange [0:0.02]
set key bottom center
set pointsize 0.5
st(x)=2*x/(1e-3*5000**3)
plot 'outputs/two_wave/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'Two Wave', 'outputs/three_wave/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'Three Wave', 'outputs/modified_three_wave/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'Modified Three Wave', 'outputs/multi_point/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'Multi Point', 'outputs/multi_point_iso/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'Multi Point Iso', 'outputs/ZB_AR1D_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB AR1D LS', 'outputs/ZB_ARMD_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB ARMD LS', 'outputs/ZB_ARMD_LSM/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB ARMD LSM', 'outputs/ZB_ARMDM_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB ARMDM LS', 'outputs/ZB_ARMDMAT_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB ARMDMAT LS', 'outputs/ZB_ARMDMMAT_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB ARMDMMAT LS', 'outputs/ZB_AM_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB AM LS', 'outputs/ZB_AM_LSM/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB AM LSM', 'outputs/ZB_AMISO_LS/coeffs_export.csv' using 7:(st($9)) smooth unique w lp lw 2 title 'ZB AMISO LS','fun3D_st.csv' using 1:2 smooth unique w l lc black lw 2 title 'LAURA','fun3D_st.csv' using (-$1):2 smooth unique w l lc black lw 2 notitle
