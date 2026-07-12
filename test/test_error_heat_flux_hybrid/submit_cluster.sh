#!/bin/bash
# Usage: bash submit_cluster.sh
# Submits one PBS job per scheme in schemes.txt

ROOT=$(pwd)

for SCHEME in $(grep -v '^#\|^$' schemes.txt); do

  JOBSCRIPT=$(mktemp /tmp/pbs_${SCHEME}_XXXX.sh)

  cat > $JOBSCRIPT << EOF
#!/bin/bash
#PBS -N CTEST_HEAT_FLUX_${SCHEME}
#PBS -o output_heat_flux_${SCHEME}.txt
#PBS -e err_heat_flux_${SCHEME}.txt
#PBS -q dicam2CPUQ
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=96:mpiprocs=96:mem=120gb
#PBS -m abe
#PBS -M vincent.delmas@hotmail.fr

source ~/.bashrc

cd \$PBS_O_WORKDIR

NPROCS=96 ctest -R "test_error_heat_flux_hybrid_${SCHEME}" 2>&1 | tee myoutput_heat_flux_${SCHEME}
EOF

  echo "Submitting job for scheme: $SCHEME"
  qsub $JOBSCRIPT
  rm $JOBSCRIPT

done
