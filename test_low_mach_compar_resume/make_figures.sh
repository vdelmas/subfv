#!/bin/bash
# Regenerate every figure locally (pvbatch cannot render on Curta: no OSMesa).
ROOT=$(cd "$(dirname "$0")"; pwd)
CASES="gresho_quad gresho_tri convergence_gresho_quad convergence_gresho_tri sedov_tri sedov_hex half_cylinder_quad_ns half_cylinder_tri_ns"
SCHEMES="multi_point three_wave multi_point_pressure WIP WIP2_NOLM"
for c in $CASES; do
  for s in $SCHEMES; do
    [ -d "$ROOT/$c/outputs/$s" ] || continue
    ( cd "$ROOT/$c" && SCHEME=$s bash visu.sh > /tmp/visu_${c}_${s}.log 2>&1 \
        && echo "ok   $c/$s" || echo "part $c/$s" )
  done
done
