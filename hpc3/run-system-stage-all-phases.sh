#!/bin/bash

# This script accept a system name as input, where the system needs to have a folder
# in the systems/ folder containing starting configurations and topology files,
# followed by a stage number
# Example call:
# $ hpc3/run-stages-123-all-phases.sh 1AZ8 1C5T

# Helper function to exit script in case of error
function fail() {
  echo -e "ERROR: $1"
  exit 1
}

# Define variables to point where we expect the different scripts to be
# Note that this expects this script to be run from the repo root directory
RUN_SCRIPT=scripts/run_simulations.sh
[ -e $RUN_SCRIPT ] || fail "Run script $RUN_SCRIPT not found"
STRUCTURES_SCRIPT=scripts/prepare_nes_structures.sh
[ -e $STRUCTURES_SCRIPT ] || fail "Structures script $STRUCTURES_SCRIPT not found"
EQUILIBRIUM_SCRIPT=hpc3/equilibrium.sh
[ -e $EQUILIBRIUM_SCRIPT ] || fail "Run script $EQUILIBRIUM_SCRIPT not found"
NES_SCRIPT=hpc3/nes.sh
[ -e $NES_SCRIPT ] || fail "Run script $NES_SCRIPT not found"
# Define the size of the simulation swarm for NES simulations
NUM_NES=100

RUN=1

system="$1"
  echo "System $system"

# Check that system directory exists
SYSTEM_DIR=systems/$system
[ -d "$SYSTEM_DIR" ] || fail "Directory $SYSTEM_DIR not found"
mkdir -p "$SYSTEM_DIR"/slurm_output

stage="$2"
# Submit equilibrium run (minimization, equilibration and production runs)
jobidEq=$(sbatch --parsable --job-name="$system"-stage$stage-equilibrium \
  --error="$SYSTEM_DIR"/slurm_output/equilibrium-s$stage.err \
  --output="$SYSTEM_DIR"/slurm_output/equilibrium-s$stage.out \
  --export=RUN_SCRIPT=$RUN_SCRIPT,SYSTEM_DIR="$SYSTEM_DIR",STAGE=$stage,STRUCTURES_SCRIPT=$STRUCTURES_SCRIPT,NUM_NES=$NUM_NES,RUN=$RUN \
  $EQUILIBRIUM_SCRIPT)

# For relevant stages, submit array of NES runs that depend on successful completion
# of the equilibrium run (i.e., will only run once the previous run is completed)
if [ "${stage:0:1}" -ne 1 ] && [ -z "${stage:2:1}" ]; then
  sbatch --dependency=afterok:"$jobidEq" --job-name="$system"-stage$stage-nes \
    --error="$SYSTEM_DIR"/slurm_output/nes-s$stage-%a.err \
    --output="$SYSTEM_DIR"/slurm_output/nes-s$stage-%a.out \
    --export=RUN_SCRIPT=$RUN_SCRIPT,SYSTEM_DIR="$SYSTEM_DIR",STAGE=$stage,RUN=$RUN \
    --array=1-$NUM_NES \
    $NES_SCRIPT
fi
