#!/bin/bash
set -e

SCHEME=$1
NPROC=${NPROCS:-4}
SECOND=$2
METHOD=$3

if [[ "$SECOND" == ".true." ]]; then
  OUTDIR=outputs/${SCHEME}_o2_${METHOD}
else
  OUTDIR=outputs/${SCHEME}
fi

mkdir -p $OUTDIR
cd $OUTDIR

sed -e "s/SCHEME_PLACEHOLDER/$SCHEME/" \
    -e "s/SECOND_ORDER_PLACEHOLDER/$SECOND/" \
    -e "s/METHOD_PLACEHOLDER/$METHOD/" \
  ../../input.template > input_data.f

gmsh -3 ../../sphere_radial.geo -o sphere_radial.msh
if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 sphere_radial.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC ../../../../build/subfvns input_data.f > log.txt

date > run.timestamp
