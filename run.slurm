#!/bin/bash
#SBATCH --job-name=waternes
#SBATCH --output=slurm.out

CONFIG=$1

echo "Running with config: $CONFIG"

bash scripts/rbfe/run.sh $CONFIG
