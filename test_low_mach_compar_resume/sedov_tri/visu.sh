#!/bin/bash
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
SCHEME=${SCHEME:-WIP}
OUTDIR=outputs/${SCHEME}

mkdir -p "${ROOT}/${OUTDIR}"
cd "${ROOT}/${OUTDIR}"

pvbatch "${ROOT}/extract_sedov_profile.py"

gnuplot -e "SCHEME='${SCHEME}'; ROOT='${ROOT}'" "${ROOT}/plot_sedov.gnu"

# screenshots need OpenGL/OSMesa, unavailable on Curta compute nodes as of
# 2026-09-15 (no swrast_dri.so, no libOSMesa) -- best effort only, run last so
# a rendering crash doesn't take down the density profile plot above.
# Rerun this script locally (where rendering works) after rsync for images/*.png.
pvbatch "${ROOT}/visualize_fields.py" || echo "WARN: pvbatch image rendering failed (no OpenGL/OSMesa on this node) -- rerun visu.sh locally after rsync to get images/*.png"

date > visu.timestamp
