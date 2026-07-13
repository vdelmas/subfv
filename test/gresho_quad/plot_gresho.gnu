if (!exists("SCHEME"))  SCHEME  = "unknown"
if (!exists("ROOT"))    ROOT    = "../.."

OUTDIR = ROOT."/outputs/".SCHEME

set terminal pngcairo enhanced font "Times,14" size 800,600
set output OUTDIR."/velocity_profile.png"

set xlabel "Radius r_{xy}"
set ylabel "|V|"
unset title

gresho(r) = (r < 0.2) ? 5.0*r : ((r < 0.4) ? 2.0 - 5.0*r : 0.0)

plot OUTDIR."/velocity_profile_initial.dat" using 1:2 with points title "Initial (t=0)", \
     OUTDIR."/velocity_profile_final.dat"   using 1:2 with points title "Final (t=10^{-2})", \
     gresho(x) with lines title "Analytical" lc rgb "black"
