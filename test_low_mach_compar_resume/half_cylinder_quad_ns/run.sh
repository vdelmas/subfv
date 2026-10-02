#!/bin/bash
# Standalone hypersonic half-cylinder (Navier-Stokes, quad mesh) test for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}
NPROC=${NPROCS:-32}

SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" "${ROOT}/input.template" > input_data.f

gmsh -3 "${ROOT}/half_cylinder_quad_ns.geo" -o half_cylinder_quad_ns.msh
if [ "$NPROC" -gt 1 ]; then
  ${SUBFV_GMSH:-subfv-gmsh} -3 half_cylinder_quad_ns.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "${BIN}/subfvns" input_data.f > log.txt

date > run.timestamp
