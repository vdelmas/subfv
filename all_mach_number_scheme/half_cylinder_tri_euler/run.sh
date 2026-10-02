#!/bin/bash
# Standalone hypersonic half-cylinder (Euler, tri mesh) test for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}
NPROC=${NPROCS:-32}

OUTDIR=outputs/WIP
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

cp "${ROOT}/input.template" input_data.f

gmsh -3 "${ROOT}/half_cylinder_tri_euler.geo" -o half_cylinder_tri_euler.msh
if [ "$NPROC" -gt 1 ]; then
  ${SUBFV_GMSH:-subfv-gmsh} -3 half_cylinder_tri_euler.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "${BIN}/subfvns" input_data.f > log.txt

date > run.timestamp
