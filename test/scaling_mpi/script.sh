#!/bin/bash

#############################
# les directives Slurm vont ici:

# Your job name (displayed by the queue)
#SBATCH -J HelloWorld

# walltime (hh:mm::ss)
#SBATCH -t 00:10:00

# Specify the number of nodes(nodes=) and the number of cores per nodes(tasks-pernode=) to be used
#SBATCH -N 1
#SBATCH --ntasks-per-node=32

# change working directory
# SBATCH --chdir=.

# fin des directives PBS
#############################

# useful informations to print
echo "#############################" 
echo "User:" $USER
echo "Date:" `date`
echo "Host:" `hostname`
echo "Directory:" `pwd`
echo "SLURM_JOBID:" $SLURM_JOBID
echo "SLURM_SUBMIT_DIR:" $SLURM_SUBMIT_DIR
echo "SLURM_JOB_NODELIST:" $SLURM_JOB_NODELIST
echo "#############################" 

#############################

# What you actually want to launch
./run.sh

# all done
echo "Job finished" 
