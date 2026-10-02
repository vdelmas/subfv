#!/bin/bash
# Standalone low-Mach convergence study (Gresho, Voronoi mesh) for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}

SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" "${ROOT}/input.template" > input_data.f

cp "${ROOT}/gresho_voro_per.cgns" .

mpirun ${MPIRUN_FLAGS} -np 1 "${BIN}/subfvgreshoconv" input_data.f > log.txt

date > run.timestamp
