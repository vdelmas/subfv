#!/bin/bash
#SBATCH -J lag_cv_voro
#SBATCH -t 06:00:00
#SBATCH -n 1
#SBATCH -p imb-resources

# Lagrangian Gresho low-Mach convergence: sweeps Ma = 1e-1..1e-4 for one
# (scheme, mesh) pair, each to t_final = 1e-2, matching exactly what
# ns_gresho_convergence_main.F90 does on the Euler side so the two sets of
# numbers are comparable. The density fluctuation itself is computed
# afterwards by conv_post.py as sqrt(sum_i V_i (rho_i - 1)^2), the same
# definition as l2_rho there.

module load gcc/11.2.0 lapack mpi blas
export MPIRUN_FLAGS="--bind-to none"

SCHEME=${SCHEME:-vp_pp}
MESHFILE=${MESHFILE:-gresho_voro.cgns}
MESHTAG=${MESHTAG:-voro}
ROOT=$SLURM_SUBMIT_DIR
BIN=/scratch/vdelmas/subfv/build_all_mach_number_scheme

for MACH in 1e-1 1e-2 1e-3 1e-4; do
  OUT=$ROOT/outputs_conv/${SCHEME}_gresho_${MESHTAG}/Ma${MACH}
  rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
  sed -e "s/MESH_PLACEHOLDER/${MESHFILE}/" \
      -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
      -e "s/MACH_PLACEHOLDER/${MACH}/" \
      -e "s/t_max = 1.0/t_max = 0.01/" \
      -e "s/n_sol_vtu = 21/n_sol_vtu = 2/" "$ROOT/input.template" > input_data.f
  cp "$ROOT/${MESHFILE}" .
  mpirun $MPIRUN_FLAGS -np 1 "$BIN/subfvlagrange" input_data.f > log.txt 2>&1
  echo "Ma=${MACH} rc=$?" >> "$ROOT/outputs_conv/${SCHEME}_gresho_${MESHTAG}/summary.txt"
done
date > "$ROOT/outputs_conv/${SCHEME}_gresho_${MESHTAG}/run.timestamp"
