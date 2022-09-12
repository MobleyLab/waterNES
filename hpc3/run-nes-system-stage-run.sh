#!/bin/bash

# This script accepts system names as input, where the system needs to have a folder
# in the systems/ folder containing starting configurations and topology files
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
NES_SCRIPT=hpc3/nes_stage_run.sh
[ -e $NES_SCRIPT ] || fail "Run script $NES_SCRIPT not found"

RUN=1
system=$1
stage=$2
runNumber=$3

echo "System $system"

# Check that system directory exists
SYSTEM_DIR=systems/$system
[ -d "$SYSTEM_DIR" ] || fail "Directory $SYSTEM_DIR not found"
mkdir -p "$SYSTEM_DIR"/slurm_output

sbatch --job-name="$system"-stage"$stage"-nes-"$runNumber" \
  --error="$SYSTEM_DIR"/slurm_output/nes-s"$stage"-"$runNumber".err \
  --output="$SYSTEM_DIR"/slurm_output/nes-s"$stage"-"$runNumber".out \
  --export=RUN_SCRIPT=$RUN_SCRIPT,SYSTEM_DIR="$SYSTEM_DIR",STAGE="$stage",RUN=$RUN,RUNNUMBER="$runNumber" \
  $NES_SCRIPT
