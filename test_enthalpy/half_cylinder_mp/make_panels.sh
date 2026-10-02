#!/bin/bash
# Per-case field images with NO colour-bar margin, plus one standalone colour bar per field,
# for the LaTeX tables in ../report. Everything here is an original render: nothing is cropped
# out of a larger image afterwards.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
unset DISPLAY

FILES=$(ls outputs/*/output_-1.pvtu 2>/dev/null)
[ -z "$FILES" ] && { echo "no finished run in outputs/"; exit 1; }
read SCHL HMIN HMAX RMIN RMAX PMIN PMAX <<< "$(pvbatch ranges.py $FILES 2>/dev/null | tail -1)"
echo "shared ranges: schlieren 0..$SCHL  H $HMIN..$HMAX  rho $RMIN..$RMAX  p $PMIN..$PMAX"

mkdir -p panels
for f in $FILES; do
  tag=$(basename "$(dirname "$f")")
  pvbatch render.py "$f" "panels/$tag" "$SCHL" "$HMIN" "$HMAX" \
          "$RMIN" "$RMAX" "$PMIN" "$PMAX" panels 2>/dev/null | grep '^saved'
done

# Colour bars: one render per field, camera off the domain so only the 2D bar overlay remains.
one=$(echo $FILES | awk '{print $1}')
pvbatch render.py "$one" "panels/bar" "$SCHL" "$HMIN" "$HMAX" \
        "$RMIN" "$RMAX" "$PMIN" "$PMAX" panels bars 2>/dev/null | grep '^saved'
