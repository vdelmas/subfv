#!/bin/bash
#SBATCH -J noh_ref
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH -p imb-resources

# Noh 2D lagrangien sur maillage raffine : quad 100x100 (10000 mailles) et tri
# non structure (10066 mailles), soit 4x le cas existant a 2500 mailles.
# L'ordre 2 de `classic` est connu pour partir en NaN sur prismes triangulaires
# (cf. le commentaire de test/lagrange_sedov_compar/job.sbatch), d'ou les deux
# ordres lances sur chaque maillage plutot que l'ordre 2 seul.

module load gcc/11.2.0 lapack mpi blas
export MPIRUN_FLAGS="--bind-to none"

SCHEME=${SCHEME:-classic}
MESH=${MESH:-noh_quad}
ORDER=${ORDER:-.true.}
TAG=${TAG:-o2}

ROOT=$SLURM_SUBMIT_DIR
BIN=/scratch/vdelmas/subfv/build_all_mach_number_scheme
OUT=$ROOT/outputs/${SCHEME}_${MESH}_${TAG}

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
sed -e "s/MESH_PLACEHOLDER/${MESH}/" -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
    -e "s/ORDER_PLACEHOLDER/${ORDER}/" "$ROOT/input.template" > input_data.f
cp "$ROOT/${MESH}.msh" .

echo "User:" $USER; echo "Date:" `date`; echo "Host:" `hostname`
echo "SLURM_JOBID:" $SLURM_JOBID; echo "scheme=${SCHEME} mesh=${MESH} order=${ORDER}"

mpirun $MPIRUN_FLAGS -np 1 "$BIN/subfvlagrange" input_data.f > log.txt 2>&1
rc=$?
echo "EXIT_CODE: $rc" >> log.txt
date > run.timestamp
exit $rc
