if (!exists("SCHEME"))  SCHEME  = "unknown"
if (!exists("ROOT"))    ROOT    = "../.."

OUTDIR = ROOT."/outputs/".SCHEME
DAT = OUTDIR."/convergence.dat"

set terminal pdf enhanced dashed font "Times,14" size 4in,3in
set output OUTDIR."/images/convergence.pdf"

# Wong (2011) colorblind-safe palette (black reserved for reference lines)
set linetype 1 lc rgb "#E69F00" lw 2 pt 7  ps 0.6 dt 1
set linetype 2 lc rgb "#56B4E9" lw 2 pt 5  ps 0.6 dt 2
set linetype 3 lc rgb "#009E73" lw 2 pt 9  ps 0.6 dt 3

set xlabel "Mach"
set ylabel "{/Symbol D}{/Symbol r}_{L^2}"
set logscale xy
set grid
set key bottom right

# Reference slopes anchored at first data point
ref_val = system("awk 'NR==2{print $2}' " . DAT) + 0.0
ma_ref  = system("awk 'NR==2{print $1}' " . DAT) + 0.0

plot DAT using 1:2 with linespoints lt 1 title "Numerical", \
     ref_val * (x/ma_ref)**(1) with lines lc rgb "black" dt 2 lw 1.5 title "slope 1", \
     ref_val * (x/ma_ref)**(2) with lines lc rgb "black" dt 3 lw 1.5 title "slope 2"
