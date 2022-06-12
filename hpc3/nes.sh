#!/bin/bash
#SBATCH --account=DMOBLEY_LAB
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00

# This script expects the following variables to be passed in via Slurm's --export flag:
# RUN_SCRIPT
# SYSTEM_DIR
# STAGE
# STRUCTURES_SCRIPT
# NUM_NES

# The source command can be replaced by loading a module (if the GROMACS version you
# want to use is available) or by a pointer to any other GROAMCS installation.
# The next line is telling shellcheck not to worry about the sourced file.
# shellcheck source=/dev/null
source ~/bin/gmx2022.1/bin/GMXRC

# Run NES simulation. The simulation number is defined by the Slurm array task ID.
bash "$RUN_SCRIPT" -d "$SYSTEM_DIR"/stage"$STAGE" -t "$SYSTEM_DIR" -x gmx \
  -o "-ntmpi $SLURM_CPUS_PER_TASK" -s "$STAGE" -p NES -n "$SLURM_ARRAY_TASK_ID"
