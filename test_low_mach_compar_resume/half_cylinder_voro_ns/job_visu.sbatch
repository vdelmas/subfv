#!/bin/bash
#SBATCH -J hcvoro_visu
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH -p visu
#SBATCH --gres=gpu:1
#SBATCH --mem=16G

# Rendu des champs sur le cluster, PNG rapatries ensuite (voir memoire
# curta-remote-rendering). Le ParaView du $HOME est le binaire EGL de Kitware :
# il n'a AUCUN repli en rendu logiciel, donc il segfault sur un noeud de calcul
# ordinaire ("failed to open swrast") et sur le noeud de connexion. Seule la
# partition visu avec --gres=gpu:1 (Quadro P4000) fonctionne. Le warning
# "X server connection. DISPLAY=" est normal et sans consequence.

module load gcc/11.2.0 lapack mpi blas
PVBATCH=$HOME/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/pvbatch

ROOT=$SLURM_SUBMIT_DIR
for S in multi_point three_wave multi_point_pressure WIP WIP2_NOLM two_wave modified_three_wave; do
  D=$ROOT/outputs/$S
  [ -f "$D/output_-1.pvtu" ] || { echo "$S : pas de sortie finale, ignore"; continue; }
  ( cd "$D" && mkdir -p images && $PVBATCH "$ROOT/visualize_fields.py" > visu_render.log 2>&1 )
  echo "$S -> $(ls $D/images/*.png 2>/dev/null | wc -l) png"
done
date > $ROOT/visu_render.timestamp
