#!/bin/bash
# Standalone shear test for the WIP scheme (extracted from test/test_shear).
# Runs the 4 mesh refinements sequentially, single rank each (matches the
# original ctest default).
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
# Dedicated build tree for this study (never the shared test/build/ used by
# other in-flight work on this machine/cluster) — see all_mach_number_scheme/README.md.
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}

MESHES="mesh_shear_1 mesh_shear_2 mesh_shear_4 mesh_shear_8"

for MESH in $MESHES; do
  OUTDIR=outputs/WIP/${MESH}
  mkdir -p "${ROOT}/${OUTDIR}"
  cd "${ROOT}/${OUTDIR}"

  sed -e "s/MESH_PLACEHOLDER/$MESH/" "${ROOT}/input.template" > input_data.f

  gmsh -3 "${ROOT}/${MESH}.geo" -o ${MESH}.msh

  mpirun ${MPIRUN_FLAGS} -np 1 "${BIN}/subfvshear" > log.txt

  date > run.timestamp
  cd "${ROOT}"
done
