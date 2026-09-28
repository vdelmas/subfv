#!/bin/bash
#SBATCH -J e_quad_multi_point_enthalpy_o1
#SBATCH -t 12:00:00
#SBATCH -n 32
#SBATCH -p imb
#SBATCH --mem-per-cpu=4G
#SBATCH -o /scratch/vdelmas/subfv/test_enthalpy/slurm-quad_multi_point_enthalpy_o1-%j.out

echo "User:" $USER; echo "Date:" `date`; echo "Host:" `hostname`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST

module purge
module load gcc/11.2.0 lapack mpi blas

HERE=/scratch/vdelmas/subfv/test_enthalpy/half_cylinder_ns
OUT=$HERE/outputs/quad_multi_point_enthalpy_o1
mkdir -p $OUT && cd $OUT

MESH=half_cylinder_quad_ns.msh
gmsh -3 $HERE/half_cylinder_quad_ns.geo -o $MESH > gmsh.log 2>&1 || { echo "gmsh FAILED"; exit 1; }

# subfv-gmsh needs GLIBC_2.29 but Curta's system glibc is 2.28 -- invoke it through the
# glibc/2.41 module's own dynamic linker (libstdc++/libgcc_s still come from the gcc module).
module load glibc/2.41
LD=$GLIBC_PATH/lib/ld-linux-x86-64.so.2
LIBPATH=$GLIBC_PATH/lib:/gpfs/softs/contrib/apps/gcc/11.2.0/lib64
$LD --library-path $LIBPATH /scratch/vdelmas/subfv/tools/subfv-gmsh/subfv-gmsh \
  -3 $MESH -part 32 -part_split -part_ghosts > gmsh_part.log 2>&1 || { echo "partition FAILED"; exit 1; }
module unload glibc/2.41

sed -e "s/MESH_PLACEHOLDER/$MESH/" \
    -e "s/SCHEME_PLACEHOLDER/multi_point_enthalpy/" \
    -e "s/SECOND_ORDER_PLACEHOLDER/.false./" \
    -e "s/METHOD_PLACEHOLDER/2/" \
    -e "s/NMAXITER_PLACEHOLDER/40000/" \
    $HERE/input.template > input_data.f

mpirun -np 32 /scratch/vdelmas/subfv/build_mp/subfvns input_data.f > log.txt 2>&1
echo "SOLVER EXIT=$?"
date > run.timestamp

# Wall Cp / heat-flux export for the LAURA comparison. Pure data extraction (no Render), so
# it needs no GL context and runs here rather than on the visu partition. Non-fatal: the run
# itself is the deliverable, the CSV can always be regenerated from coeffs_-1.pvtu.
PVPYTHON=$HOME/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvpython
if [ -f coeffs_-1.pvtu ]; then
  $PVPYTHON $HERE/export_coeffs.py coeffs_-1.pvtu coeffs_export.csv > export.log 2>&1 \
    || echo "coeffs export FAILED (see export.log) -- rerun locally"
fi
tail -20 log.txt
