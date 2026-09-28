#!/bin/bash
# Render the schlieren + total-enthalpy panel for every run in outputs/, all on one shared
# colour range so the four panels (quad/tri x three_wave/three_wave_enthalpy) are directly
# comparable. DISPLAY must be unset: pvbatch renders fully offscreen.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
unset DISPLAY

FILES=$(ls outputs/*/output_-1.pvtu 2>/dev/null)
[ -z "$FILES" ] && { echo "no finished run in outputs/"; exit 1; }

read SCHL HMIN HMAX <<< "$(pvbatch ranges.py $FILES 2>/dev/null | tail -1)"
echo "shared ranges: schlieren 0..$SCHL   H $HMIN..$HMAX"

mkdir -p figures
for f in $FILES; do
  tag=$(basename "$(dirname "$f")")
  pvbatch render.py "$f" "figures/$tag" "$SCHL" "$HMIN" "$HMAX" 2>/dev/null | grep '^saved'
done
