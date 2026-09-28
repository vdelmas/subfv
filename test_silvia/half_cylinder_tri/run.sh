#!/bin/bash
set -e

ORDER=${1:-1}
NPROC=${NPROCS:-1}

OUTDIR=outputs/o${ORDER}

mkdir -p $OUTDIR
cd $OUTDIR

# génération input
sed -e "s/ORDER_PLACEHOLDER/$ORDER/" \
  ../../input.template > input_data.f

# exécution
gmsh -3 ../../half_cylinder_tri.geo -o half_cylinder_tri.msh
if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 half_cylinder_tri.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC ../../../../build/subfvsilvia input_data.f > log.txt

date > run.timestamp
