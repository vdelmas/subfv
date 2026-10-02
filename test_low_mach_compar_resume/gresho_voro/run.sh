#!/bin/bash
# Standalone Gresho vortex (Voronoi mesh) test for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}
NPROC=1  # polyhedral CGNS + periodic stitching: serial only

SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" "${ROOT}/input.template" > input_data.f

# periodic mesh requires single proc (same limitation as test/gresho_quad)
if grep -q "periodic_mesh *= *\.true\." input_data.f; then
  NPROC=1
fi

cp "${ROOT}/gresho_voro_per.cgns" .
mpirun ${MPIRUN_FLAGS} -np $NPROC "${BIN}/subfvns" input_data.f > log.txt

date > run.timestamp
