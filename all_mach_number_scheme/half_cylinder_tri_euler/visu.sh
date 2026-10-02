#!/bin/bash
set -e

OUTDIR=outputs/WIP
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# screenshots need OpenGL/OSMesa, unavailable on Curta compute nodes as of
# 2026-09-15 (no swrast_dri.so, no libOSMesa) -- best effort only. Rerun this
# script locally (where rendering works) after rsync to get images/*.png.
pvbatch ../../visualize_fields.py || echo "WARN: pvbatch image rendering failed (no OpenGL/OSMesa on this node) -- rerun visu.sh locally after rsync to get images/*.png"

date > visu.timestamp
