#!/bin/bash
# Submits the WIP2 sweep (low-Mach-scaled nodal pressure) to Slurm on Curta, for
# the 8 cases relevant to the 4 WIP2 acceptance criteria (Gresho velocity
# preservation, density-fluctuation-vs-Mach convergence, Sedov robustness,
# Mach-27 half-cylinder wall heat flux). Does NOT touch the existing WIP
# baseline results (outputs/WIP/) already present from the previous campaign.
#
# Run this FROM the cluster copy of this directory
# (e.g. /scratch/vdelmas/subfv/all_mach_number_scheme/), after building
# ../build_all_mach_number_scheme (see README.md).
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)

CASES="gresho_quad gresho_tri convergence_gresho_quad convergence_gresho_tri sedov_tri sedov_hex half_cylinder_quad_ns half_cylinder_tri_ns"
SCHEMES="multi_point_pressure"

for CASE in $CASES; do
  for SCHEME in $SCHEMES; do
    echo "== submitting $CASE / $SCHEME =="
    # regenerate this scheme's own output subdir from scratch every time
    # (never trust a stale partitioned mesh, see curta-cluster-config memory)
    # -- outputs/WIP and any other scheme's outputs are left untouched.
    rm -rf "${ROOT}/${CASE}/outputs/${SCHEME}"
    ( cd "${ROOT}/${CASE}" && sbatch -J "zb_${CASE}_${SCHEME}" --export=ALL,SCHEME="${SCHEME}" job.sbatch )
  done
done
