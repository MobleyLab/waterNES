#!/bin/bash

function fail()
{
  echo -e "ERROR: $1"
  exit 1
}

# We expect this to run from the root of the repository
[ -e tests/$(basename $0) ] || fail "Expected to be in the root repository."

RUN_SCRIPT=scripts/run_simulations.sh
# Check that we can find the script we want to test
[ -e $RUN_SCRIPT ] || fail "Can't find script $RUN_SCRIPT."

# Check that GROMACS is available as `gmx`
GMX=gmx
which $GMX &> /dev/null || fail "Executable $GMX not found."

# Do some changes to input files to ensure fast runs
for file in input/*.mdp; do cp $file $file.bk; done
perl -pi -e 's/^nsteps *=\K.+/ 2/' input/*.mdp

TEST_SYSTEM=1HPX
WORKDIR=systems/$TEST_SYSTEM
TOPDIR=$WORKDIR
CRD=$WORKDIR/system.gro
# Check that typical usage runs without errors, using very small simulations to be fast
bash $RUN_SCRIPT -d $WORKDIR/sys2 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 2 -p min || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/sys2 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 2 -p eqNVT || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/sys2 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 2 -p eqNPT || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/sys2 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 2 -p prod || fail "Error in simulation."
# Run the other systems in one go
bash $RUN_SCRIPT -d $WORKDIR/sys1 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 1 -p "min eqNVT eqNPT prod" || fail "Error in simulation."
bash $RUN_SCRIPT -d $WORKDIR/sys3 -t $TOPDIR -c $CRD -x $GMX -o "-nt 1" -s 3 -p "min eqNVT eqNPT prod" || fail "Error in simulation."

# Restore the original mdp files
for file in input/*.mdp; do mv $file.bk $file; done

# If we reached here, all commands passed without error!
echo "SUCCESS: All simulations ran without error."
exit 0
