#!/bin/bash
set -e

SCHEME=$1
NPROC=${NPROCS:-4}
SECOND=${2:-.false.}
METHOD=${3:-0}

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

gmsh -3 ../../forward_facing_quad.geo -o forward_facing_quad.msh
if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 forward_facing_quad.msh -part $NPROC -part_split -part_ghosts
fi

mpirun ${MPIRUN_FLAGS} -np $NPROC ../../../../build/subfvns input_data.f > log.txt

date > run.timestamp
