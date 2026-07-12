#!/bin/bash
# Génère les 4 cas aux coins pour la visualisation du maillage.
# 1 seul CPU, 10 itérations, schéma multi_point_iso.
# Sorties dans corners_visu/{label}/output_0.pvtu

set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
GEO=half_cylinder_tri_ns_error_heat
CORNERS_DIR=${ROOT}/corners_visu
SCHEME=multi_point_iso

mkdir -p ${CORNERS_DIR}

declare -A CORNERS
CORNERS["top_left"]="2.0 2.0000e-05"      # fsize=2.0  blsize=2e-5   (FF grossier, CL fine)
CORNERS["top_right"]="2.0 1.0000e-03"     # fsize=2.0  blsize=1e-3   (FF grossier, CL grossière)
CORNERS["bottom_left"]="0.2 2.0000e-05"   # fsize=0.2  blsize=2e-5   (FF fin,     CL fine)
CORNERS["bottom_right"]="0.2 1.0000e-03"  # fsize=0.2  blsize=1e-3   (FF fin,     CL grossière)

for LABEL in top_left top_right bottom_left bottom_right; do
  read FSIZE BLSIZE <<< "${CORNERS[$LABEL]}"
  OUTDIR=${CORNERS_DIR}/${LABEL}
  MESHDIR=${ROOT}/meshes/${FSIZE}_${BLSIZE}

  echo "=== ${LABEL}: fsize=${FSIZE}, blsize=${BLSIZE} ==="

  # ---- Maillage ----
  if [ ! -f "${MESHDIR}/${GEO}.msh" ]; then
    mkdir -p ${MESHDIR}
    sed -e "s/^fsize=.*;/fsize=${FSIZE};/" \
        -e "s/^blsize=.*;/blsize=${BLSIZE};/" \
      ${ROOT}/${GEO}.geo > ${MESHDIR}/${GEO}.geo
    gmsh -3 ${MESHDIR}/${GEO}.geo -o ${MESHDIR}/${GEO}.msh -format msh4 2>/dev/null
    echo "  mesh generated: ${MESHDIR}/${GEO}.msh"
  else
    echo "  mesh already exists"
  fi

  # ---- Input ----
  mkdir -p ${OUTDIR}
  sed -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
      -e "s/SECOND_ORDER_PLACEHOLDER/.false./" \
      -e "s/METHOD_PLACEHOLDER/1/" \
      -e "s|meshfile_path='',|meshfile_path='${MESHDIR}/',|" \
      -e "s/half_cylinder_tri_ns\.msh/${GEO}.msh/" \
      -e "s/n_max_iter=80000/n_max_iter=10/" \
      -e "s/n_iter_print=8000/n_iter_print=1/" \
      -e "s/n_iter_write_sol=999999/n_iter_write_sol=999999/" \
    ${ROOT}/input.template > ${OUTDIR}/input_data.f

  # ---- Solver ----
  cd ${OUTDIR}
  mpirun -np 1 ${ROOT}/../../build/subfvns input_data.f > log.txt 2>&1
  echo "  solver done -> $(ls output_0.pvtu 2>/dev/null && echo OK || echo MISSING)"
  cd ${ROOT}

done

echo ""
echo "Done. Files in: ${CORNERS_DIR}"
echo "Run pvbatch visu_corners.py to generate mesh images."
