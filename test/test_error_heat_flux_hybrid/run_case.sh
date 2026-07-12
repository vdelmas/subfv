#!/bin/bash
set -e

SCHEME=$1
FSIZE=$2
BLSIZE=$3
NPROC=${NPROCS:-1}
SECOND_ORDER=${SECOND_ORDER:-.false.}
METHOD=${METHOD:-1}

ROOT=$(cd "$(dirname "$0")"; pwd)
GEO=half_cylinder_tri_ns_error_heat
MESHDIR=${ROOT}/meshes/${FSIZE}_${BLSIZE}

# Generate mesh only if not already done
if [ ! -f "${MESHDIR}/${GEO}.msh" ]; then
  mkdir -p ${MESHDIR}
  sed -e "s/^fsize=.*;/fsize=${FSIZE};/" \
      -e "s/^blsize=.*;/blsize=${BLSIZE};/" \
    ${ROOT}/${GEO}.geo > ${MESHDIR}/${GEO}.geo
  gmsh -3 ${MESHDIR}/${GEO}.geo -o ${MESHDIR}/${GEO}.msh -format msh4 2>/dev/null
fi

# Partition mesh only if not already done for this proc count
if [ "$NPROC" -gt 1 ] && [ ! -f "${MESHDIR}/${GEO}_${NPROC}/${GEO}_1.msh" ]; then
  mkdir -p ${MESHDIR}/${GEO}_${NPROC}
  subfv-gmsh -3 ${MESHDIR}/${GEO}.msh -part $NPROC -part_split -part_ghosts -format msh4 \
    -o ${MESHDIR}/${GEO}_${NPROC}/${GEO}.msh 2>/dev/null
fi

mkdir -p ${ROOT}/outputs/${SCHEME}/${FSIZE}_${BLSIZE}
cd ${ROOT}/outputs/${SCHEME}/${FSIZE}_${BLSIZE}

if [ "$NPROC" -gt 1 ]; then
  MESHFILE_PATH="${MESHDIR}/${GEO}_${NPROC}/"
else
  MESHFILE_PATH="${MESHDIR}/"
fi

sed -e "s/SCHEME_PLACEHOLDER/$SCHEME/" \
    -e "s/SECOND_ORDER_PLACEHOLDER/$SECOND_ORDER/" \
    -e "s/METHOD_PLACEHOLDER/$METHOD/" \
    -e "s|meshfile_path='',|meshfile_path='${MESHFILE_PATH}',|" \
    -e "s/half_cylinder_tri_ns\.msh/${GEO}.msh/" \
  ${ROOT}/input.template > input_data.f

mpirun ${MPIRUN_FLAGS} -np $NPROC ${ROOT}/../../build/subfvns input_data.f > log.txt

date > run.timestamp
