#!/bin/bash
# Build one summary PDF per scheme from the consolidated campaign in
# ../test_low_mach_compar_resume/ (all schemes run with ONE binary, so the
# numbers in the five reports are directly comparable).
#
# The substitution is done in python, not sed: sed reinterprets the backslash
# in a "\_" replacement and silently drops it, which put the scheme name into
# math mode and mangled the whole title line.
ROOT=$(cd "$(dirname "$0")"; pwd)
SCHEMES="${@:-multi_point three_wave multi_point_pressure WIP WIP2_NOLM}"

for S in $SCHEMES; do
  python3 - "$ROOT" "$S" <<'PY'
import sys
root, s = sys.argv[1], sys.argv[2]
tpl = open(f"{root}/report_template.tex").read()
out = tpl.replace("@SCHEME@", s.replace("_", r"\_")).replace("@SCHEMEPATH@", s)
open(f"{root}/report_{s}.tex", "w").write(out)
PY
  ( cd "$ROOT" && pdflatex -interaction=nonstopmode "report_${S}.tex" > "build_${S}.log" 2>&1
                  pdflatex -interaction=nonstopmode "report_${S}.tex" > "build_${S}.log" 2>&1 )
  P=$(pdfinfo "$ROOT/report_${S}.pdf" 2>/dev/null | awk '/Pages/{print $2}')
  M=$(pdftotext "$ROOT/report_${S}.pdf" - 2>/dev/null | grep -c 'non disponible')
  printf '%-34s %s pages, %s encadre(s) vide(s)\n' "report_${S}.pdf" "$P" "$M"
done
