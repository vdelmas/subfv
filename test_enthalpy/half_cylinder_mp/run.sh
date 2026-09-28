#!/bin/bash
# Inviscid M=20 half cylinder run through subfvns, for the multidimensional (node-based)
# solvers: multi_point (= the paper's Gallice-2D) and multi_point_enthalpy (= MGallice-2D).
# Same meshes and same freestream as ../half_cylinder, so the panels are directly comparable
# to the three_wave / three_wave_enthalpy ones.
#
# usage: ./run.sh <quad|tri> <multi_point|multi_point_enthalpy|three_wave|three_wave_enthalpy>
set -e

MESHKIND=$1
SCHEME=$2
NPROC=${NPROCS:-6}
# Enthalpy preservation is an exact property of the CONVERGED steady state only: the discrete
# steady state with h == h_inf satisfies the energy equation automatically once mass is
# satisfied, but nothing forces h = h_inf during the transient. So the measured departure keeps
# shrinking with the residual and the iteration count is what sets the quality of the result --
# at 12000/5000 iterations (the ~5 min-per-case calibration, ~41 iter/s on quad and ~17 iter/s
# on tri at 6 ranks with two cases side by side) the quad was still an order of magnitude off
# machine precision and the tri had not even finished establishing (||d rho||/||rho|| = 0.19).
# These defaults target a genuinely settled solution instead; always confirm with steadiness.py.
# Calibrated for ~10 min per case at 6 ranks with ONE case running at a time: measured
# 109 iter/s on the quad mesh alone, and the tri mesh has 2.35x the cells. Running cases
# concurrently is a false economy here -- 4 at once on this box drops the quad to 30 iter/s,
# i.e. less than half the aggregate throughput of running them one after another.
if [ -z "$NMAXITER" ]; then
  if [ "$MESHKIND" = "quad" ]; then NMAXITER=65000; else NMAXITER=28000; fi
fi

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$HERE/../..
OUTDIR=$HERE/outputs/${MESHKIND}_${SCHEME}

rm -rf "$OUTDIR"; mkdir -p "$OUTDIR"; cd "$OUTDIR"

cp "$HERE/cyl_${MESHKIND}.msh" .
sed -e "s/MESH_PLACEHOLDER/cyl_${MESHKIND}.msh/" \
    -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
    -e "s/NMAXITER_PLACEHOLDER/${NMAXITER}/" \
    "$HERE/input.template" > input_data.f

if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 cyl_${MESHKIND}.msh -part $NPROC -part_split -part_ghosts > gmsh_part.log 2>&1
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC "$ROOT/build/subfvns" input_data.f > log.txt 2>&1

date > run.timestamp
