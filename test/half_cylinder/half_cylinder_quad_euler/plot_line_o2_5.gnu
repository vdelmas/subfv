set terminal pdf enhanced font 'Times,12' size 4in,3in
set linetype 1 lc rgb "#E69F00" lw 2 pt 7
set linetype 2 lc rgb "#56B4E9" lw 2 pt 9
set linetype 3 lc rgb "#009E73" lw 2 pt 5
set linetype 4 lc rgb "#F0E442" lw 2 pt 11
set linetype 5 lc rgb "#0072B2" lw 2 pt 13
set linetype 6 lc rgb "#D55E00" lw 2 pt 6
set linetype 7 lc rgb "#CC79A7" lw 2 pt 8
set linetype 8 lc rgb "#000000" lw 2 pt 4
set output 'plot_line_o2_5.pdf'
set datafile separator ','
set grid
set xlabel 'x'
set ylabel 'Density'
set xrange [-2:-1]
set key top left
set pointsize 0.5
plot 
