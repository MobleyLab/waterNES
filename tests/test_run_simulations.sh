#!/bin/bash

function restore_input_files()
{
  # Restore the original mdp files
  for file in input/*.mdp; do [ -e "$file".bk ] && mv "$file".bk "$file"; done
}

function fail()
{
  restore_input_files
  echo -e "ERROR: $1"
  exit 1
}

# We expect this to run from the root of the repository
[ -e tests/"$(basename "$0")" ] || fail "Expected to be in the root repository."

RUN_SCRIPT=scripts/run_simulations.sh
STRUCTURES_SCRIPT=scripts/prepare_nes_structures.sh
# Check that we can find the script we want to test
[ -e $RUN_SCRIPT ] || fail "Can't find script $RUN_SCRIPT."
[ -e $STRUCTURES_SCRIPT ] || fail "Can't find script $STRUCTURES_SCRIPT."

# Check that GROMACS is available as `gmx`
GMX=gmx
which $GMX &> /dev/null || fail "Executable $GMX not found."

# Do some changes to input files to ensure fast runs
for file in input/*.mdp; do cp "$file" "$file".bk; done
perl -pi -e 's/^nsteps *=\K.+/ 4/' input/*.mdp
perl -pi -e 's/^nstxout *=\K.+/ 2/' input/*.mdp

TEST_SYSTEM=1HPX
WORKDIR=systems/$TEST_SYSTEM
TOPDIR=$WORKDIR
CRD=$WORKDIR/system.gro
# Check that typical usage runs without errors, using very small simulations to be fast
bash $RUN_SCRIPT -d $WORKDIR/stage2 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 2 -p min || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/stage2 -t $TOPDIR -x $GMX -o "-nt 1" -s 2 -p eqNVT || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/stage2 -t $TOPDIR -x $GMX -o "-nt 1" -s 2 -p eqNPT || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/stage2 -t $TOPDIR -x $GMX -o "-nt 1" -s 2 -p prod || fail "Error in simulation."
# Run the other stages in one go
bash $RUN_SCRIPT -d $WORKDIR/stage1 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 1 -p "min eqNVT eqNPT prod" || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/stage3 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 3 -p "min eqNVT eqNPT prod" || fail "Error in simulation."

# Prepare NES structures for stages 2 and 3, using different numbers of structures
bash $STRUCTURES_SCRIPT -d $WORKDIR/stage2 -x $GMX -n 2 || fail "Error in NES structure preparation"
bash $STRUCTURES_SCRIPT -d $WORKDIR/stage3 -x $GMX -n 1 || fail "Error in NES structure preparation"

# Run the NES stages
bash $RUN_SCRIPT -d $WORKDIR/stage2 -t $TOPDIR -x $GMX -o "-nt 1" -s 2 -p NES -n 1 || fail "Error in NES simulation."
bash $RUN_SCRIPT -d $WORKDIR/stage2 -t $TOPDIR -x $GMX -o "-nt 1" -s 2 -p NES -n 2 || fail "Error in NES simulation."
bash $RUN_SCRIPT -d $WORKDIR/stage3 -t $TOPDIR -x $GMX -o "-nt 1" -s 3 -p NES -n 1 || fail "Error in NES simulation."

# If we reached here, all commands passed without error!
restore_input_files
echo "SUCCESS: All simulations ran without error."
exit 0
