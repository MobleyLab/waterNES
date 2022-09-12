#!/bin/bash

for system in "$@"; do
  SYSTEM_DIR=systems/$system
  sbatch --job-name="$system"-distances \
      --error="$SYSTEM_DIR"/slurm_output/distances.err \
      --output="$SYSTEM_DIR"/slurm_output/distances.out \
      --export=SYSTEM="$system" \
      hpc3/compute-distances.sh
done
