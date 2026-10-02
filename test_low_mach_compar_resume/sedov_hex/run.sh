#!/bin/bash
# Standalone spherical Sedov blast test (hex mesh) for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}
NPROC=${NPROCS:-32}

SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" "${ROOT}/input.template" > input_data.f

gmsh -3 "${ROOT}/spherical_sedov_hex.geo" -o spherical_sedov_hex.msh
if [ "$NPROC" -gt 1 ]; then
  ${SUBFV_GMSH:-subfv-gmsh} -3 spherical_sedov_hex.msh -part $NPROC -part_split -part_ghosts
fi

mpirun ${MPIRUN_FLAGS} -np $NPROC "${BIN}/subfvns" input_data.f > log.txt

date > run.timestamp
