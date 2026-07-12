#!/bin/bash
set -e

: "${NPROCS:=4}"

if ! [[ "$NPROCS" =~ ^[0-9]+$ ]] || [ "$NPROCS" -lt 1 ]; then
    echo "Erreur: NPROCS invalide: $NPROCS" >&2
    exit 1
fi

NPS=$(for ((p=1; p<=NPROCS; p*=2)); do echo -n "$p "; done)
echo $NPS

if [ -f scaling_mpi.txt ]; then
    rm scaling_mpi.txt
fi

gmsh -3 1drp.geo -o 1drp.msh
for NP in $NPS;do
  subfv-gmsh -3 1drp.msh -part $NP -part_split -part_ghosts
  mpirun ${MPIRUN_FLAGS} -np $NP ../../build/subfvns input_data.f > log_${NP}.txt
  cpu_time=$(awk '/CPU Time \(s\)/ {print $4}' log_${NP}.txt)
  echo "$NP $cpu_time" >> scaling_mpi.txt
done

date > run.timestamp
