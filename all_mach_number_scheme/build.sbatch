#!/bin/bash
#SBATCH -J amn_build
#SBATCH -t 00:20:00
#SBATCH -n 4
#SBATCH -p imb

module load gcc/11.2.0 lapack mpi blas

cd /scratch/vdelmas/subfv
mkdir -p build_all_mach_number_scheme
cd build_all_mach_number_scheme
# Les trois _ROOT sont indispensables dans un repertoire de build NEUF : sans eux
# cmake ne trouve ni CGNS ni HDF5 ni ParMETIS et construit SILENCIEUSEMENT sans le
# lecteur .cgns parallele (binaire ~421 Ko au lieu de ~1.1 Mo, maillages de Voronoi
# illisibles). Un build existant marche sans eux parce que son CMakeCache les a deja.
# Controle : cmake.log doit afficher "CGNS reader with ParMETIS partitioning".
cmake -DCMAKE_BUILD_TYPE=Release \
  -DCGNS_ROOT=/gpfs/home/vdelmas/opt/cgns \
  -DPARMETIS_ROOT=/gpfs/home/vdelmas/opt/parmetis \
  -DHDF5_ROOT=/gpfs/softs/contrib/apps/hdf5/1.14.4/openmpi_5.0.2/gcc_11.2.0 \
  .. > cmake.log 2>&1
make -j4 subfvns subfvshear subfvgreshoconv subfvlagrange > make.log 2>&1
rc=$?
echo "BUILD_EXIT_CODE: $rc"
date > build.timestamp
exit $rc
