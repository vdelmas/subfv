#!/bin/bash
# Post-processing for the standalone shear test (polar error plot per mesh).
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)

MESHES="mesh_shear_1 mesh_shear_2 mesh_shear_4 mesh_shear_8"

for MESH in $MESHES; do
  OUTDIR=outputs/WIP/${MESH}
  cd "${ROOT}/${OUTDIR}"
  gnuplot -e "SCHEME='WIP'; MESH='${MESH}'" "${ROOT}/plot_shear.gnu"
  date > visu.timestamp
  cd "${ROOT}"
done
