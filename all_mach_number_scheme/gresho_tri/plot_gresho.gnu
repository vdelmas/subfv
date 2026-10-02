if (!exists("SCHEME"))  SCHEME  = "unknown"
if (!exists("ROOT"))    ROOT    = "../.."

OUTDIR = ROOT."/outputs/".SCHEME

set terminal pdf enhanced dashed font "Times,14" size 4in,3in
set output OUTDIR."/images/velocity_profile.pdf"

# Wong (2011) colorblind-safe palette (black reserved for analytical solutions)
set linetype 1 lc rgb "#E69F00" lw 2 pt 7  ps 0.6 dt 1
set linetype 2 lc rgb "#56B4E9" lw 2 pt 5  ps 0.6 dt 2
set linetype 3 lc rgb "#009E73" lw 2 pt 9  ps 0.6 dt 3
set linetype 4 lc rgb "#0072B2" lw 2 pt 13 ps 0.6 dt 4
set linetype 5 lc rgb "#D55E00" lw 2 pt 11 ps 0.6 dt 5
set linetype 6 lc rgb "#CC79A7" lw 2 pt 6  ps 0.6 dt 6

set xlabel "Radius r_{xy}"
set ylabel "|V|"
unset title

gresho(r) = (r < 0.2) ? 5.0*r : ((r < 0.4) ? 2.0 - 5.0*r : 0.0)

plot OUTDIR."/velocity_profile_final.dat" using 1:2 with points lt 1 title "Numerical solution", \
     gresho(x) with lines lc rgb "black" title "Analytical"
