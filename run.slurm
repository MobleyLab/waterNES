#!/bin/bash
#SBATCH --job-name=waternes
#SBATCH --output=slurm.out

CONFIG=$1

echo "Running with config: $CONFIG"

bash src/waternes/pipelines/rbfe/run_simulations.sh $CONFIG
