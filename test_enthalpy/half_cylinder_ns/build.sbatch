#!/bin/bash
#SBATCH -J ent_build
#SBATCH -t 00:30:00
#SBATCH -n 8
#SBATCH -p imb
#SBATCH -o /scratch/vdelmas/subfv/test_enthalpy/slurm-build-%j.out

echo "User:" $USER; echo "Date:" `date`; echo "Host:" `hostname`
module purge
module load gcc/11.2.0 lapack mpi blas

# Fresh out-of-source dir: /scratch/vdelmas/subfv/build holds a CMakeCache.txt that was
# rsynced from the local workstation and points at the local source path, so cmake refuses
# to reuse it. Build into build_ent instead of clobbering that directory.
cd /scratch/vdelmas/subfv
rm -rf build_ent
mkdir -p build_ent && cd build_ent
cmake /scratch/vdelmas/subfv -DCMAKE_BUILD_TYPE=Release
make -j8 subfvns subfveulerho
echo "=== BUILD RESULT ==="
ls -la /scratch/vdelmas/subfv/build_ent/subfvns /scratch/vdelmas/subfv/build_ent/subfveulerho
