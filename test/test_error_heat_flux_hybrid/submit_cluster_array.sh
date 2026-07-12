#!/bin/bash
# Usage: bash submit_cluster_array.sh [SCHEME]
# Submits one PBS job array per scheme, one array element per fsize x blsize case

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

  JOBSCRIPT=$(mktemp /tmp/pbs_array_${SCHEME}_XXXX.sh)

  # PBS array header: use -J only if more than 1 case
  if [ $NCASES -gt 1 ]; then
    ARRAY_DIRECTIVE="#PBS -J 0-${N}"
    INDEX_VAR="\$PBS_ARRAY_INDEX"
  else
    ARRAY_DIRECTIVE=""
    INDEX_VAR="0"
  fi

  cat > $JOBSCRIPT << EOF
#!/bin/bash
#PBS -N HEAT_FLUX_${SCHEME}
#PBS -o output_heat_flux_${SCHEME}_^array_index^.txt
#PBS -e err_heat_flux_${SCHEME}_^array_index^.txt
#PBS -q dicam2CPUQ
#PBS -l walltime=08:00:00
#PBS -l select=4:ncpus=24:mpiprocs=24:mem=30gb
${ARRAY_DIRECTIVE}
#PBS -m abe
#PBS -M vincent.delmas@hotmail.fr

source ~/.bashrc

cd \$PBS_O_WORKDIR

CASES=($(printf '"%s" ' "${CASES[@]}"))

read FSIZE BLSIZE <<< "\${CASES[${INDEX_VAR}]}"

MPIRUN_FLAGS="-machinefile \$PBS_NODEFILE" NPROCS=96 bash run_case.sh ${SCHEME} \$FSIZE \$BLSIZE 2>&1 | tee log_array_${INDEX_VAR}.txt
EOF

  echo "Submitting ${NCASES} case(s) for scheme: $SCHEME"
  qsub $JOBSCRIPT
  rm $JOBSCRIPT

done
