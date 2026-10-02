#!/bin/bash
# Consolidated low-Mach / hypersonic comparison for the subfvns schemes.
#
# Runs the whole acceptance suite -- Gresho vortex, low-Mach density-fluctuation
# convergence, Sedov robustness, Mach-27 half-cylinder wall heat flux -- for
# every scheme we compare, so that all the numbers quoted anywhere come from one
# single campaign with one single binary.
#
# Run from the cluster copy of this directory, after building
# ../build_all_mach_number_scheme (see ../all_mach_number_scheme/README.md).
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)

CASES="gresho_quad gresho_tri gresho_voro convergence_gresho_quad convergence_gresho_tri convergence_gresho_voro sedov_tri sedov_hex half_cylinder_quad_ns half_cylinder_tri_ns half_cylinder_voro_ns"
SCHEMES="${SCHEMES:-multi_point three_wave multi_point_pressure WIP WIP2_NOLM}"

for CASE in $CASES; do
  for SCHEME in $SCHEMES; do
    echo "== $CASE / $SCHEME =="
    # always regenerate mesh+partition from the .geo (see curta-cluster-config:
    # a stale partition file is indistinguishable from a real MPI bug)
    rm -rf "${ROOT}/${CASE}/outputs/${SCHEME}"
    ( cd "${ROOT}/${CASE}" && sbatch -J "lm_${CASE}_${SCHEME}" \
        --export=ALL,SCHEME="${SCHEME}" job.sbatch )
  done
done
