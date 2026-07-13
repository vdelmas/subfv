if (!exists("SCHEME"))  SCHEME  = "unknown"
if (!exists("ROOT"))    ROOT    = "../.."

OUTDIR = ROOT."/outputs/".SCHEME

set terminal pngcairo enhanced font "Times,14" size 800,600
set output OUTDIR."/density_profile.png"

set xlabel "Radius r_{xy}"
set ylabel "Density"
set xrange [0:1.2]
unset title
set key top left inside offset 1, 0

plot ROOT."/analytic_sedov_2D.dat" using 2:3 with lines lw 2 lc rgb "black" title "Analytic", \
     OUTDIR."/density_profile.dat"  using 1:2 with points pt 7 ps 0.3 lc rgb "#E69F00" title SCHEME noenhanced
