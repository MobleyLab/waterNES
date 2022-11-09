#!/bin/bash

for system in "$@"; do
  SYSTEM_DIR=systems/$system
  sbatch --job-name="$system"-analysis \
      --error="$SYSTEM_DIR"/slurm_output/analysis.err \
      --output="$SYSTEM_DIR"/slurm_output/analysis.out \
      --export=SYSTEM="$system" \
      hpc3/do-postprocessing.sh
done
