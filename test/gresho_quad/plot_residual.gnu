if (!exists("OUTDIR")) OUTDIR = "."

set terminal pdf enhanced dashed font "Times,14" size 4in,3in
set output OUTDIR."/images/residual.pdf"

# Wong (2011) colorblind-safe palette (black reserved for analytical solutions)
set linetype 1 lc rgb "#E69F00" lw 2 pt 7  ps 0.6 dt 1
set linetype 2 lc rgb "#56B4E9" lw 2 pt 5  ps 0.6 dt 2
set linetype 3 lc rgb "#009E73" lw 2 pt 9  ps 0.6 dt 3
set linetype 4 lc rgb "#0072B2" lw 2 pt 13 ps 0.6 dt 4
set linetype 5 lc rgb "#D55E00" lw 2 pt 11 ps 0.6 dt 5
set linetype 6 lc rgb "#CC79A7" lw 2 pt 6  ps 0.6 dt 6

set xlabel "Iteration"
set ylabel "Residual / R_0"
set logscale y
set grid
set key top right

rho_0   = system("awk 'NR==1{print $4}' " . OUTDIR . "/residual.dat") + 0.0
rhoUx_0 = system("awk 'NR==1{print $5}' " . OUTDIR . "/residual.dat") + 0.0
rhoUy_0 = system("awk 'NR==1{print $6}' " . OUTDIR . "/residual.dat") + 0.0
rhoE_0  = system("awk 'NR==1{print $8}' " . OUTDIR . "/residual.dat") + 0.0

nlines = system("wc -l < " . OUTDIR . "/residual.dat") + 0
if (nlines > 1) {
  DAT = "< head -n -1 ".OUTDIR."/residual.dat"
} else {
  DAT = OUTDIR."/residual.dat"
}

plot DAT using 1:($4/rho_0)   with linespoints lt 1 title "{/Symbol r}", \
     DAT using 1:($5/rhoUx_0) with linespoints lt 2 title "{/Symbol r}u_x", \
     DAT using 1:($6/rhoUy_0) with linespoints lt 3 title "{/Symbol r}u_y", \
     DAT using 1:($8/rhoE_0)  with linespoints lt 5 title "{/Symbol r}E"
