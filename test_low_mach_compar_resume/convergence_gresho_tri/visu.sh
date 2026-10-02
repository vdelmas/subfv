#!/bin/bash
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}

mkdir -p "${ROOT}/${OUTDIR}/images"
cd "${ROOT}/${OUTDIR}"

gnuplot -e "SCHEME='${SCHEME}'; ROOT='${ROOT}'" "${ROOT}/plot_convergence.gnu"

date > visu.timestamp
