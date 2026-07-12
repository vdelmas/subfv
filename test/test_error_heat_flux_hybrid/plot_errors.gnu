if (!exists("SCHEME")) SCHEME = "multi_point_iso"

set terminal pdfcairo enhanced font "Times,56" size 16in, 13in
set datafile separator ' '
set palette defined (0 "#440154", 1 "#31688e", 2 "#35b779", 3 "#fde725")
set size ratio 1
set xlabel "First cell size" font "Times,60"
set ylabel "Far field cell size" font "Times,60"
set yrange [0.1:2.1]
set xtics ("2×10^{-5}" -4.69897, \
           "3×10^{-5}" -4.51020, \
           "5×10^{-5}" -4.32142, \
           "7×10^{-5}" -4.13264, \
           "10^{-4}"   -3.94386, \
           "2×10^{-4}" -3.75510, \
           "3×10^{-4}" -3.56633, \
           "4×10^{-4}" -3.37755, \
           "6×10^{-4}" -3.18877, \
           "10^{-3}"   -3.00000) \
    font "Times,56" rotate by -45
set ytics font "Times,50"
set cbtics font "Times,50"

set format cb "10^{%.0f}"
if (exists("RES_MIN") && exists("RES_MAX")) set cbrange [log10(RES_MIN):log10(RES_MAX)]
set output sprintf("residuals_%s.pdf", SCHEME)
set title sprintf("Final residual -- %s", SCHEME) noenhanced font "Times,56"
set cblabel "R_{final} / R_0" font "Times,56"
plot sprintf("residuals_%s.dat", SCHEME) \
  using (log10($2)):1:(log10($3)) with image pixels notitle
unset output
set cbrange [*:*]

set format cb "%g"
if (exists("CH_MIN") && exists("CH_MAX")) set cbrange [CH_MIN:CH_MAX]
set output sprintf("errors_ch_%s.pdf", SCHEME)
set title sprintf("Error -- Ch (%s)", SCHEME) noenhanced font "Times,56"
set cblabel "Relative error Ch (%)" font "Times,56"
plot sprintf("errors_%s.dat", SCHEME) \
  using (log10($2)):1:3 with image pixels notitle
unset output
set cbrange [*:*]

if (exists("CP_MIN") && exists("CP_MAX")) set cbrange [CP_MIN:CP_MAX]
set output sprintf("errors_cp_%s.pdf", SCHEME)
set title sprintf("Error -- Cp (%s)", SCHEME) noenhanced font "Times,56"
set cblabel "Relative error Cp (%)" font "Times,56"
plot sprintf("errors_%s.dat", SCHEME) \
  using (log10($2)):1:4 with image pixels notitle
unset output
set cbrange [*:*]
