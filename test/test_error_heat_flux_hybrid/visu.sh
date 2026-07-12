#!/bin/bash
set -e

SCHEME=$1
ROOT=$(pwd)

for FSIZE in $(grep -v '^#\|^$' fsize.txt); do
  for BLSIZE in $(grep -v '^#\|^$' blsize.txt); do
    CASENAME=${FSIZE}_${BLSIZE}
    OUTDIR=outputs/${SCHEME}/${CASENAME}

    pushd $OUTDIR > /dev/null

    pvbatch ${ROOT}/extract_coeffs.py 2>/dev/null
    python3 ${ROOT}/compute_error.py > error.txt
    gnuplot -e "ROOTDIR='${ROOT}'; LABEL='${CASENAME}'" ${ROOT}/plot_case.gnu

    echo "$(cat error.txt)  ${CASENAME}"
    popd > /dev/null
  done
done
