#!/bin/bash
set -e

SCHEME=$1

ROOT=$(cd "$(dirname "$0")"; pwd)
OUTDIR=outputs/${SCHEME}

mkdir -p "$OUTDIR"
cd "$OUTDIR"

sed -e "s/SCHEME_PLACEHOLDER/$SCHEME/" \
  ../../input.template > input_data.f

gmsh -3 ../../gresho_quad.geo -o gresho_quad.msh

mpirun ${MPIRUN_FLAGS} -np 1 ../../../../build/subfvgreshoconv input_data.f > log.txt

date > run.timestamp
