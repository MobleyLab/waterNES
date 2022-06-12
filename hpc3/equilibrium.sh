#!/bin/bash
#SBATCH --account=DMOBLEY_LAB
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00

# This script expects the following variables to be passed in via Slurm's --export flag:
# RUN_SCRIPT
# SYSTEM_DIR
# STAGE

# The source command can be replaced by loading a module (if the GROMACS version you
# want to use is available) or by a pointer to any other GROMACS installation.
# The next line is telling shellcheck not to worry about the sourced file.
# shellcheck source=/dev/null
source ~/bin/gmx2022.1/bin/GMXRC

# Run equilibrium simulations (minimization, equilibration, production)
bash "$RUN_SCRIPT" -d "$SYSTEM_DIR"/stage"$STAGE" -t "$SYSTEM_DIR" \
  -c "$SYSTEM_DIR"/system.gro -x gmx -o "-ntmpi $SLURM_CPUS_PER_TASK" \
  -s "$STAGE" -p "min eqNVT eqNPT prod" || exit 1

# For stages 2 and 3, prepare the structures for NES simulations after completion
# of the production run. This allows to start the NES simulations as a dependency.
if [ "$STAGE" -eq 2 ] || [ "$STAGE" -eq 3 ]; then
  bash "$STRUCTURES_SCRIPT" -d "$SYSTEM_DIR"/stage"$STAGE" -x gmx -n "$NUM_NES"
fi
