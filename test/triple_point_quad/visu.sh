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

pvbatch "${ROOT}/visualize_fields.py"

date > visu.timestamp
