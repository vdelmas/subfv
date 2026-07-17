#!/bin/bash
set -e

SCHEME=$1

ROOT=$(cd "$(dirname "$0")"; pwd)
OUTDIR=outputs/${SCHEME}

mkdir -p "${ROOT}/${OUTDIR}/images"
cd "${ROOT}/${OUTDIR}"

gnuplot -e "SCHEME='${SCHEME}'; ROOT='${ROOT}'" "${ROOT}/plot_convergence.gnu"

date > visu.timestamp
