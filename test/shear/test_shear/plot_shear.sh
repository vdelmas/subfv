#!/bin/bash
set -e

SCHEME=$1
MESH=$2

cd outputs/${SCHEME}/${MESH}
gnuplot -e "SCHEME='${SCHEME}'; MESH='${MESH}'" ../../../plot_shear.gnu
