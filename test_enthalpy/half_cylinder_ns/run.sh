#!/bin/bash
# Hypersonic M_inf = 17.6 half cylinder, VISCOUS (subfvns), isothermal wall at 500 K.
#
# Same geometry/freestream as the project's existing test/half_cylinder_{quad,tri}_ns cases;
# only the Riemann solver changes, so three_wave and three_wave_enthalpy are directly
# comparable. h_inf = gamma/(gamma-1) p/rho + |u|^2/2 = 3.5*57615.576818 + 0.5*5000^2
#         = 1.2701654e7.
# Unlike the Euler case, total enthalpy is NOT expected to be uniform here: viscous work and
# wall heat flux make h drop inside the boundary layer. That is the point -- with
# three_wave_enthalpy the h deficit that remains is physical (boundary layer) rather than the
# numerical undershoot three_wave leaves across the bow shock, which is exactly what a
# total-enthalpy-based boundary-layer criterion needs.
#
# usage: ./run.sh <quad|tri> <three_wave|three_wave_enthalpy> [second_order] [method]
set -e

MESHKIND=$1
SCHEME=$2
SECOND=${3:-.false.}
METHOD=${4:-2}
NPROC=${NPROCS:-4}

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$HERE/../..
if [[ "$SECOND" == ".true." ]]; then
  TAG=${MESHKIND}_${SCHEME}_o2_m${METHOD}
else
  TAG=${MESHKIND}_${SCHEME}_o1
fi
OUTDIR=$HERE/outputs/$TAG

rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"; cd "$OUTDIR"

MESH=half_cylinder_${MESHKIND}_ns.msh
gmsh -3 "$HERE/half_cylinder_${MESHKIND}_ns.geo" -o "$MESH" > gmsh.log 2>&1
sed -e "s/MESH_PLACEHOLDER/${MESH}/" \
    -e "s/NMAXITER_PLACEHOLDER/${NMAXITER:-40000}/" \
    -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
    -e "s/SECOND_ORDER_PLACEHOLDER/${SECOND}/" \
    -e "s/METHOD_PLACEHOLDER/${METHOD}/" \
    "$HERE/input.template" > input_data.f

if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 "$MESH" -part $NPROC -part_split -part_ghosts > gmsh_part.log 2>&1
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "$ROOT/build/subfvns" input_data.f > log.txt 2>&1

date > run.timestamp
