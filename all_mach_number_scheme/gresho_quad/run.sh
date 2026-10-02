#!/bin/bash
# Standalone Gresho vortex (quad mesh) test for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}
NPROC=${NPROCS:-4}

SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" "${ROOT}/input.template" > input_data.f

# periodic mesh requires single proc (same limitation as test/gresho_quad)
if grep -q "periodic_mesh *= *\.true\." input_data.f; then
  NPROC=1
fi

gmsh -3 "${ROOT}/gresho_quad.geo" -o gresho_quad.msh
if [ "$NPROC" -gt 1 ]; then
  ${SUBFV_GMSH:-subfv-gmsh} -3 gresho_quad.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "${BIN}/subfvns" input_data.f > log.txt

date > run.timestamp
