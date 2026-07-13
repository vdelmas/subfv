#!/bin/bash
set -e

SCHEME=$1
SECOND=${2:-.false.}
METHOD=${3:-0}

if [[ "$SECOND" == ".true." ]]; then
  OUTDIR=outputs/${SCHEME}_o2_${METHOD}
else
  OUTDIR=outputs/${SCHEME}
fi

ROOT=$(cd "$(dirname "$0")"; pwd)

mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

pvbatch "${ROOT}/extract_sedov_profile.py"

pvbatch "${ROOT}/visualize_fields.py"

gnuplot -e "SCHEME='${SCHEME}'; ROOT='${ROOT}'" "${ROOT}/plot_sedov.gnu"

date > visu.timestamp
