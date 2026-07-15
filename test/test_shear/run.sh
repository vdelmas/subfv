#!/bin/bash
set -e

SCHEME=$1
MESH=$2
NPROC=${NPROCS:-1}

OUTDIR=outputs/${SCHEME}/${MESH}

mkdir -p $OUTDIR
cd $OUTDIR

sed -e "s/SCHEME_PLACEHOLDER/$SCHEME/" \
    -e "s/MESH_PLACEHOLDER/$MESH/" \
  ../../../input.template > input_data.f

gmsh -3 ../../../${MESH}.geo -o ${MESH}.msh
if [ "$NPROC" -gt 1 ]; then
  subfv-gmsh -3 ${MESH}.msh -part $NPROC -part_split -part_ghosts
fi

mpirun ${MPIRUN_FLAGS} -np $NPROC ../../../../../build/subfvshear > log.txt

date > run.timestamp
