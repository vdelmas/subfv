#!/bin/bash
#SBATCH -J lag_gr_voro
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH -p imb-resources

# Lagrangian Gresho vortex, order 1, on a VORONOI mesh (polyhedral CGNS read
# block-wise and partitioned in-code by src/partition/, no subfv-gmsh split).
# Walltime is twice the quad/tri job's: the Voronoi cells have a smaller
# inradius, so the CFL-limited dt is ~1.9x smaller for the same cell count.

module load gcc/11.2.0 lapack mpi blas
export MPIRUN_FLAGS="--bind-to none"

SCHEME=${SCHEME:-classic}
MESHFILE=${MESHFILE:-gresho_voro.cgns}
MESHTAG=${MESHTAG:-voro}
MACH=${MACH:-1e-2}
TMAX=${TMAX:-0.25}

ROOT=$SLURM_SUBMIT_DIR
BIN=/scratch/vdelmas/subfv/build_all_mach_number_scheme
OUT=$ROOT/outputs/${SCHEME}_gresho_${MESHTAG}

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"

sed -e "s/MESH_PLACEHOLDER/${MESHFILE}/" \
    -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
    -e "s/MACH_PLACEHOLDER/${MACH}/" \
    -e "s/t_max = 1.0/t_max = ${TMAX}/" "$ROOT/input.template" > input_data.f
cp "$ROOT/${MESHFILE}" .

echo "User:" $USER ; echo "Date:" `date` ; echo "Host:" `hostname`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "scheme=${SCHEME} mesh=${MESHFILE} mach=${MACH} tmax=${TMAX}"

mpirun $MPIRUN_FLAGS -np 1 "$BIN/subfvlagrange" input_data.f > log.txt 2>&1
rc=$?
echo "EXIT_CODE: $rc" >> log.txt
date > run.timestamp
exit $rc
