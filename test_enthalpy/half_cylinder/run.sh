#!/bin/bash
# Hypersonic M=20 half cylinder, aho Euler solver (subfveulerho).
# usage: ./run.sh <quad|tri> <three_wave|three_wave_enthalpy> [order] [tmax]
set -e

MESHKIND=$1
FLUX=$2
ORDER=${3:-1}
TMAX=${4:-8.0}
NPROC=${NPROCS:-4}

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$HERE/../..
OUTDIR=$HERE/outputs/${MESHKIND}_${FLUX}_o${ORDER}

rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"; cd "$OUTDIR"

cp "$HERE/cyl_${MESHKIND}.msh" .
sed -e "s/MESH_PLACEHOLDER/cyl_${MESHKIND}.msh/" \
    -e "s/FLUX_PLACEHOLDER/${FLUX}/" \
    -e "s/ORDER_PLACEHOLDER/${ORDER}/" \
    -e "s/TMAX_PLACEHOLDER/${TMAX}/" \
    "$HERE/input.template" > input_data.f

if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 cyl_${MESHKIND}.msh -part $NPROC -part_split -part_ghosts > gmsh_part.log 2>&1
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "$ROOT/build/subfveulerho" input_data.f > log.txt 2>&1

date > run.timestamp
