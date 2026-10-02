#!/bin/bash
# Standalone Gresho vortex (tri mesh) test for the WIP scheme.
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
BIN=${SUBFV_BUILD:-${ROOT}/../../build_all_mach_number_scheme}
NPROC=${NPROCS:-32}

SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}
mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" "${ROOT}/input.template" > input_data.f

# periodic mesh requires single proc: periodic-boundary stitching
# (mesh_connectivity_module.F90's periodic stencil build) requires matching
# left/right boundary face counts on EACH rank's local partition, which
# arbitrary MPI partitioning does not guarantee -- confirmed on Curta at
# NPROC=32 ("Error periodic boundaries !" / error stop). Same limitation
# gresho_quad/run.sh already guards against; this was a real gap here.
if grep -q "periodic_mesh *= *\.true\." input_data.f; then
  NPROC=1
fi

gmsh -3 "${ROOT}/gresho_tri.geo" -o gresho_tri.msh
if [ "$NPROC" -gt 1 ]; then
  ${SUBFV_GMSH:-subfv-gmsh} -3 gresho_tri.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "${BIN}/subfvns" input_data.f > log.txt

date > run.timestamp
