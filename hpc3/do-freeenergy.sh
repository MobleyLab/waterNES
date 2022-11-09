#!/bin/bash
#SBATCH --account=DMOBLEY_LAB
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00

# This script expects the following variables to be passed in via Slurm's --export flag:
# SYSTEM

source /data/homezvol2/pmerz/.bash_profile
module load anaconda
conda activate waterNES

python scripts/analyze.py "$SYSTEM" systems/"$SYSTEM"/run1 --read --doFreeEnergyOfStages 1 4 7 6 6.1 6.2 6.3 6.4 6.5 6.6 6.7 5
