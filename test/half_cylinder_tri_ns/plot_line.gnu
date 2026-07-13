set terminal pdf enhanced font 'Times,12' size 4in,3in
set linetype 1 lc rgb "#E69F00" lw 2 pt 7
set linetype 2 lc rgb "#56B4E9" lw 2 pt 9
set linetype 3 lc rgb "#009E73" lw 2 pt 5
set linetype 4 lc rgb "#F0E442" lw 2 pt 11
set linetype 5 lc rgb "#0072B2" lw 2 pt 13
set linetype 6 lc rgb "#D55E00" lw 2 pt 6
set linetype 7 lc rgb "#CC79A7" lw 2 pt 8
set linetype 8 lc rgb "#000000" lw 2 pt 4
set output 'plot_line.pdf'
set datafile separator ','
set grid
set xlabel 'x'
set ylabel 'Density'
set xrange [-2:-1]
set key top left
set pointsize 0.5
plot 'outputs/two_wave/line_export.csv' using 19:22 w lp lw 2 title 'Two Wave', 'outputs/three_wave/line_export.csv' using 19:22 w lp lw 2 title 'Three Wave', 'outputs/modified_three_wave/line_export.csv' using 19:22 w lp lw 2 title 'Modified Three Wave', 'outputs/multi_point/line_export.csv' using 19:22 w lp lw 2 title 'Multi Point', 'outputs/multi_point_iso/line_export.csv' using 19:22 w lp lw 2 title 'Multi Point Iso', 'outputs/ZB_AR1D_LS/line_export.csv' using 19:22 w lp lw 2 title 'ZB AR1D LS', 'outputs/ZB_AR1D_LSM/line_export.csv' using 19:22 w lp lw 2 title 'ZB AR1D LSM', 'outputs/ZB_ARMDM_LS/line_export.csv' using 19:22 w lp lw 2 title 'ZB ARMDM LS', 'outputs/ZB_ARMDM_LSM/line_export.csv' using 19:22 w lp lw 2 title 'ZB ARMDM LSM', 'outputs/ZB_ARMDMAT_LS/line_export.csv' using 19:22 w lp lw 2 title 'ZB ARMDMAT LS', 'outputs/ZB_AM_LS/line_export.csv' using 19:22 w lp lw 2 title 'ZB AM LS', 'outputs/ZB_AM_LSM/line_export.csv' using 19:22 w lp lw 2 title 'ZB AM LSM', 'outputs/ZB_AMISO_LS/line_export.csv' using 19:22 w lp lw 2 title 'ZB AMISO LS', 'outputs/ZB_AMISO_LSM/line_export.csv' using 19:22 w lp lw 2 title 'ZB AMISO LSM'
