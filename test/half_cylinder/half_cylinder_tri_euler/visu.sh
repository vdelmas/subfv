#!/bin/bash
set -e

SCHEME=$1
NPROC=$2
SECOND=$3
METHOD=$4

if [[ "$SECOND" == ".true." ]]; then
  OUTDIR=outputs/${SCHEME}_o2_${METHOD}
else
  OUTDIR=outputs/${SCHEME}
fi

mkdir -p $OUTDIR
cd $OUTDIR

pvbatch ../../visualize_fields.py

date > visu.timestamp
