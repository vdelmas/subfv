#!/bin/bash
#SBATCH -J lsv_voro
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH -p imb-resources

# Sedov lagrangien sur maillage de VORONOI (CGNS polyedrique, lu par bloc et
# partitionne dans le code). Walltime genereux : le dt limite par CFL est plus
# petit sur cellules polygonales (rayon inscrit plus faible), et le Sedov tri
# de reference a deja demande 3h51.

module load gcc/11.2.0 lapack mpi blas
export MPIRUN_FLAGS="--bind-to none"

SCHEME=${SCHEME:-classic}
MESHFILE=${MESHFILE:-sedov_voro.cgns}
MESHTAG=${MESHTAG:-voro}

ROOT=$SLURM_SUBMIT_DIR
BIN=/scratch/vdelmas/subfv/build_all_mach_number_scheme
OUT=$ROOT/outputs/${SCHEME}_sedov_2d_${MESHTAG}

rm -rf "$OUT"; mkdir -p "$OUT"; cd "$OUT"
sed -e "s/MESH_PLACEHOLDER/${MESHFILE}/" -e "s/SCHEME_PLACEHOLDER/${SCHEME}/" \
    "$ROOT/input.template" > input_data.f
cp "$ROOT/${MESHFILE}" .

echo "User:" $USER; echo "Date:" `date`; echo "Host:" `hostname`
echo "SLURM_JOBID:" $SLURM_JOBID; echo "scheme=${SCHEME} mesh=${MESHFILE}"

mpirun $MPIRUN_FLAGS -np 1 "$BIN/subfvlagrange" input_data.f > log.txt 2>&1
rc=$?
echo "EXIT_CODE: $rc" >> log.txt
date > run.timestamp
exit $rc
