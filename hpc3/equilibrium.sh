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
# RUN
# STRUCTURES_SCRIPT
# NUM_NES

# The source command can be replaced by loading a module (if the GROMACS version you
# want to use is available) or by a pointer to any other GROMACS installation.
# The next line is telling shellcheck not to worry about the sourced file.
# shellcheck source=/dev/null
source ~/bin/gmx2022.2+/bin/GMXRC

# Run equilibrium simulations (minimization, equilibration, production)
bash "$RUN_SCRIPT" -d "$SYSTEM_DIR"/run"$RUN"/stage"$STAGE" -t "$SYSTEM_DIR" \
  -c "$SYSTEM_DIR"/system.gro -x gmx -o "-ntomp $SLURM_CPUS_PER_TASK" \
  -s "$STAGE" -p "min eqNVT eqNPT prod" -m 3000000 || exit 1

if [ "${STAGE:0:1}" -ne 1 ] && [ -z "${STAGE:2:1}" ]; then
    bash "$STRUCTURES_SCRIPT" -d "$SYSTEM_DIR"/run"$RUN"/stage"$STAGE" -x gmx -n "$NUM_NES" -m 3000000
fi
