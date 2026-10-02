#!/bin/bash
#SBATCH -J amn_gresho_tri
#SBATCH -t 00:45:00
#SBATCH -n 1
#SBATCH -p imb

module load gcc/11.2.0 lapack mpi blas
module load gnuplot/5.4.1
export PATH="$HOME/gmsh-4.15.2-Linux64/bin:$HOME/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin:$HOME/.local/bin:$PATH"

# subfv-gmsh is a prebuilt binary needing GLIBC 2.29+, newer than Curta's
# system glibc (2.28) — invoke it through the glibc/2.41 module's own
# ld-linux, see curta-cluster-config memory.
module load glibc/2.41
GMSH_LD="$GLIBC_PATH/lib/ld-linux-x86-64.so.2"
GMSH_LIBPATH="$GLIBC_PATH/lib:/gpfs/softs/contrib/apps/gcc/11.2.0/lib64"
SUBFVGMSH_BIN="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)/tools/subfv-gmsh/subfv-gmsh"
export SUBFV_GMSH="$GMSH_LD --library-path $GMSH_LIBPATH $SUBFVGMSH_BIN"

ulimit -s unlimited

echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID

bash run.sh
rc=$?
if [ $rc -eq 0 ]; then
  bash visu.sh
  rc=$?
fi
exit $rc
