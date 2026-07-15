#!/bin/bash
set -e

SCHEME=$1

OUTDIR=outputs/periodic_${SCHEME}
mkdir -p $OUTDIR
cd $OUTDIR

# periodic BC requires single proc
gmsh -3 ../../gresho_quad.geo -o gresho_quad.msh
mpirun ${MPIRUN_FLAGS} -np 1 ../../../../build/subfvns input_data.f > log.txt

date > run.timestamp
