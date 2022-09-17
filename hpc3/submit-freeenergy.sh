#!/bin/bash

for system in "$@"; do
  SYSTEM_DIR=systems/$system
  sbatch --job-name="$system"-freeenergy \
      --error="$SYSTEM_DIR"/slurm_output/freeenergy.err \
      --output="$SYSTEM_DIR"/slurm_output/freeenergy.out \
      --export=SYSTEM="$system" \
      hpc3/do-freeenergy.sh
done
