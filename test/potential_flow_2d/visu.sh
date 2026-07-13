#!/bin/bash
set -e

SCHEME=$1
SECOND=$2
METHOD=$3

if [[ "$SECOND" == ".true." ]]; then
  OUTDIR=outputs/${SCHEME}_o2_${METHOD}
else
  OUTDIR=outputs/${SCHEME}
fi

cd $OUTDIR
mkdir -p images
pvbatch ../../visualize_fields.py
gnuplot -e "OUTDIR='.'" ../../plot_residual.gnu
