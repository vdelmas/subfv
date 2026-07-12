if (!exists("SCHEME")) SCHEME = "unknown"
if (!exists("MESH"))   MESH   = "unknown"

set terminal pdfcairo enhanced font "Times,12" size 4in,4in
set output sprintf("error_shear_%s_%s.pdf", SCHEME, MESH)

##############################
# Supprimer axes cartésiens
##############################
unset border
unset xtics
unset ytics
set zeroaxis lc rgb "white"

##############################
# Calcul du rayon max
##############################
unset polar
stats "error_shear.dat" using 2 nooutput name 'D'
max_radius = D_max * 1.1
if (max_radius <= 0) { max_radius = 1 }

step = max_radius / 5.0

# Format compact sans zéro leading dans l'exposant (1e-5 au lieu de 1e-05)
cexp(x) = (substr(sprintf("%.0e",x), strlen(sprintf("%.0e",x))-1, strlen(sprintf("%.0e",x))-1) eq "0" ? \
           substr(sprintf("%.0e",x),1,strlen(sprintf("%.0e",x))-2) . \
           substr(sprintf("%.0e",x),strlen(sprintf("%.0e",x)),strlen(sprintf("%.0e",x))) : \
           sprintf("%.0e",x))

set rtics (cexp(step) step, cexp(2*step) 2*step, cexp(3*step) 3*step, \
           cexp(4*step) 4*step, cexp(5*step) 5*step) font "Times,9"

##############################
# Mode polaire
##############################
set polar
set angles degrees
set size square
set grid polar lw 0.1 lc rgb "gray50" lt 2
set rrange [0:max_radius]
set trange [0:360]
set style data linespoints

set raxis

##############################
# Labels 0, 90, 180, 270
##############################
r_label = max_radius * 1.1
set xrange [-1.15*max_radius : 1.15*max_radius]
set yrange [-1.15*max_radius : 1.15*max_radius]

set label "0°"   at first  r_label,       0  center front
set label "90°"  at first        0, r_label  center front
set label "180°" at first -r_label,       0  center front
set label "270°" at first       0, -r_label  center front

set pointsize 0.25

##############################
# Plot
##############################
plot "error_shear.dat" using 1:2 with linespoints pt 7 lw 1.5 notitle

unset output
