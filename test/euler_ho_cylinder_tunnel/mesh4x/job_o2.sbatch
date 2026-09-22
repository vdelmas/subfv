#!/bin/bash
#SBATCH -J c4x_o2
#SBATCH -t 00:45:00
#SBATCH -n 32
#SBATCH -p imb
#SBATCH --exclude=n054,n040

module purge
module load gcc/11.2.0 lapack mpi blas
BIN=/scratch/vdelmas/subfv/build_shock_bubble_maxlambda/subfveulerho
mpirun -np 32 --bind-to none $BIN input_data_o2.f > run.log 2>&1
rc=$?
echo "EXIT_CODE: $rc" >> run.log
date > run.timestamp
exit $rc
