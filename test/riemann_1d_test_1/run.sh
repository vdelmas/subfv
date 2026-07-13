#!/bin/bash
set -e

SCHEME=$1
NPROC=${NPROCS:-4}

OUTDIR=outputs/${SCHEME}

mkdir -p $OUTDIR
cd $OUTDIR

# génération input
sed -e "s/SCHEME_PLACEHOLDER/$SCHEME/" \
  ../../input.template > input_data.f

gmsh -3 ../../1drp.geo -o 1drp.msh
if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 1drp.msh -part $NPROC -part_split -part_ghosts
fi
mpirun ${MPIRUN_FLAGS} -np $NPROC ../../../../../build/subfvns input_data.f > log.txt

date > run.timestamp
