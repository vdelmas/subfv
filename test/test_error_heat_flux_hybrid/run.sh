#!/bin/bash
set -e

SCHEME=$1
NPROC=${NPROCS:-1}
SECOND_ORDER=${SECOND_ORDER:-.false.}
METHOD=${METHOD:-1}

ROOT=$(pwd)
GEO=half_cylinder_tri_ns_error_heat

for FSIZE in $(grep -v '^#\|^$' fsize.txt); do
  for BLSIZE in $(grep -v '^#\|^$' blsize.txt); do
    CASENAME=${FSIZE}_${BLSIZE}
    OUTDIR=outputs/${SCHEME}/${CASENAME}
    mkdir -p $OUTDIR

    pushd $OUTDIR > /dev/null

    sed -e "s/^fsize=.*;/fsize=${FSIZE};/" \
        -e "s/^blsize=.*;/blsize=${BLSIZE};/" \
      ${ROOT}/${GEO}.geo > ${GEO}.geo

    gmsh -3 ${GEO}.geo -o ${GEO}.msh -format msh2 2>/dev/null

    if [ "$NPROC" -gt 1 ]; then
      subfv-gmsh -3 ${GEO}.msh -part $NPROC -part_split -part_ghosts 2>/dev/null
    fi

    sed -e "s/SCHEME_PLACEHOLDER/$SCHEME/" \
        -e "s/SECOND_ORDER_PLACEHOLDER/$SECOND_ORDER/" \
        -e "s/METHOD_PLACEHOLDER/$METHOD/" \
        -e "s/half_cylinder_tri_ns\.msh/${GEO}.msh/" \
      ${ROOT}/input.template > input_data.f

    mpirun ${MPIRUN_FLAGS} -np $NPROC ../../../../../build/subfvns input_data.f > log.txt

    date > run.timestamp
    popd > /dev/null
  done
done
