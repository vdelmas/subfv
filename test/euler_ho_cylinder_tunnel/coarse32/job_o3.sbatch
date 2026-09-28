#!/bin/bash
#SBATCH -J cyl_prod_o3c
#SBATCH -t 01:30:00
#SBATCH -n 32
#SBATCH -p imb
#SBATCH --exclude=n054,n040

echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`

module purge
module load gcc/11.2.0 lapack mpi blas

BIN=/scratch/vdelmas/subfv/build_shock_bubble_maxlambda/subfveulerho

t0=$(date +%s)
mpirun -np 32 --bind-to none $BIN input_data_o3.f > run.log 2>&1
rc=$?
t1=$(date +%s)
echo "ELAPSED_SECONDS: $((t1-t0))" >> run.log
echo "EXIT_CODE: $rc" >> run.log
date > run.timestamp
exit $rc
