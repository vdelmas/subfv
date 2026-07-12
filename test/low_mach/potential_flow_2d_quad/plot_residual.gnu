if (!exists("OUTDIR")) OUTDIR = "."

set terminal pngcairo enhanced font "Times,14" size 800,600
set output OUTDIR."/images/residual.png"

set linetype 1 lc rgb "#E69F00" lw 2
set linetype 2 lc rgb "#56B4E9" lw 2
set linetype 3 lc rgb "#009E73" lw 2
set linetype 4 lc rgb "#0072B2" lw 2

set xlabel "Iteration"
set ylabel "Residual"
set logscale y
set grid
set key top right

plot OUTDIR."/residual.dat" using 1:4 with lines title "{/Symbol r}", \
     OUTDIR."/residual.dat" using 1:5 with lines title "{/Symbol r}u_x", \
     OUTDIR."/residual.dat" using 1:6 with lines title "{/Symbol r}u_y", \
     OUTDIR."/residual.dat" using 1:8 with lines title "{/Symbol r}E"
