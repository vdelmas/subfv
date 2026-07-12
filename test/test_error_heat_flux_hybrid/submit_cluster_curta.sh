#!/bin/bash
# Usage: bash submit_cluster_curta.sh [SCHEME]
# Submits one Slurm job per scheme in schemes.txt

ROOT=$(pwd)
SCHEMES=${1:-$(grep -v '^#\|^$' schemes.txt)}

for SCHEME in $SCHEMES; do

  JOBSCRIPT=$(mktemp /tmp/slurm_${SCHEME}_XXXX.sh)

  cat > $JOBSCRIPT << EOF
#!/bin/bash
#SBATCH -J HEAT_FLUX_${SCHEME}
#SBATCH -o output_heat_flux_${SCHEME}.txt
#SBATCH -e err_heat_flux_${SCHEME}.txt
#SBATCH -p imb
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=96
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vincent.delmas@hotmail.fr

source ~/.bashrc

cd \$SLURM_SUBMIT_DIR

NPROCS=96 ctest -R "test_error_heat_flux_hybrid_${SCHEME}" 2>&1 | tee myoutput_heat_flux_${SCHEME}
EOF

  echo "Submitting job for scheme: $SCHEME"
  sbatch $JOBSCRIPT
  rm $JOBSCRIPT

done
