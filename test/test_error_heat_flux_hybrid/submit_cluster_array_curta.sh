#!/bin/bash
# Usage: bash submit_cluster_array_curta.sh [SCHEME]
# Submits one Slurm job array per scheme, one array element per fsize x blsize case

ROOT=$(pwd)
SCHEMES=${1:-$(grep -v '^#\|^$' schemes.txt)}

for SCHEME in $SCHEMES; do

  # Build ordered list of all fsize_blsize cases
  CASES=()
  while IFS= read -r FSIZE; do
    [[ -z "$FSIZE" || "$FSIZE" == \#* ]] && continue
    while IFS= read -r BLSIZE; do
      [[ -z "$BLSIZE" || "$BLSIZE" == \#* ]] && continue
      CASES+=("${FSIZE} ${BLSIZE}")
    done < blsize.txt
  done < fsize.txt

  N=$(( ${#CASES[@]} - 1 ))
  NCASES=${#CASES[@]}

  JOBSCRIPT=$(mktemp /tmp/slurm_array_${SCHEME}_XXXX.sh)

  cat > $JOBSCRIPT << EOF
#!/bin/bash
#SBATCH -J HEAT_FLUX_${SCHEME}
#SBATCH -o log_array_%a.txt
#SBATCH -e err_array_%a.txt
#SBATCH -p imb,imb-resources
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=40G
#SBATCH --array=0-${N}
#SBATCH --requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vincent.delmas@hotmail.fr

source ~/.bashrc

cd \$SLURM_SUBMIT_DIR

CASES=($(printf '"%s" ' "${CASES[@]}"))

read FSIZE BLSIZE <<< "\${CASES[\$SLURM_ARRAY_TASK_ID]}"

NPROCS=32 bash run_case.sh ${SCHEME} \$FSIZE \$BLSIZE 2>&1 | tee log_array_\${SLURM_ARRAY_TASK_ID}.txt
EOF

  echo "Submitting ${NCASES} case(s) for scheme: $SCHEME"
  sbatch $JOBSCRIPT
  rm $JOBSCRIPT

done
