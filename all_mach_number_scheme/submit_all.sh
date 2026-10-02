#!/bin/bash
# Submits all all_mach_number_scheme cases to Slurm on Curta.
# Run this FROM the cluster copy of this directory
# (e.g. /scratch/vdelmas/subfv/all_mach_number_scheme/), after building
# ../build_all_mach_number_scheme (see README.md).
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)

CASES="test_shear gresho_quad gresho_tri convergence_gresho_quad convergence_gresho_tri sedov_tri sedov_hex half_cylinder_quad_euler half_cylinder_quad_ns half_cylinder_tri_euler half_cylinder_tri_ns"

for CASE in $CASES; do
  echo "== submitting $CASE =="
  # regenerate mesh/output dir from scratch every time (never trust a stale
  # partitioned mesh, see curta-cluster-config memory)
  rm -rf "${ROOT}/${CASE}/outputs"
  ( cd "${ROOT}/${CASE}" && sbatch job.sbatch )
done

echo
echo "gresho_summary is NOT submitted here: it needs pdflatex, which Curta"
echo "does not have. Run all_mach_number_scheme/gresho_summary/combine.sh"
echo "locally after rsync'ing the 4 gresho/convergence_gresho results back."
