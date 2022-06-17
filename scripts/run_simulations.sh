#!/bin/bash

function usage()
{
  echo "USAGE: $(basename "$0") [-h] -d dir -t dir [-c crd] -x gmx [-o opts] -s N -p {min,eqNVT,eqNPT,prod,NES} [-n N]"
}

function help()
{
  echo "Run Water NES simulations."
  echo
  usage
  echo "Options:"
  echo -e "-h\tPrint this help and exit"
  echo -e "-d dir\tThe base working directory for the stage"
  echo -e "-t top\tThe topology directory including topology files"
  echo -e "-c crd\tThe input coordinates"
  echo -e "\tMandatory for the minimization phase, ignored otherwise"
  echo -e "-x gmx\tThe GROMACS executable to use"
  echo -e "-o opts\tThe run options, e.g. \"-nt 16 -ntomp 4\", optional"
  echo -e "-s N\tThe stage to run, one of:"
  echo -e "\t1: Stage with fully interacting, non-restrained water"
  echo -e "\t2: Stage with fully interacting, restrained water"
  echo -e "\t3: Stage with non-interacting, restrained water"
  echo -e "-p {min,eqNVT,eqNPT,prod,NES}"
  echo -e "\tThe simulation phase to run, one or more of: min, eqNVT, eqNPT, prod, NES"
  echo -e "\tSeparate different phases by spaces, e.g. -p \"min eqNVT eqNPT\""
  echo -e "\tNote: The NES phase is only available for stages 2 and 3"
  echo -e "-n N\tThe run number for the NES phase"
  echo -e "\tMandatory for the NES phase, ignored otherwise"
  echo
}

function fail()
{
  echo -e "ERROR: $1"
  exit 1
}

# Get the options
while getopts "d:t:c:x:o:s:p:n:h" option; do
  case $option in
    d) BASEDIR="${OPTARG}";;
    t) TOPDIR="${OPTARG}";;
    c) CRD="${OPTARG}";;
    x) GMX="${OPTARG}";;
    o) RUN_PARAMS="${OPTARG}";;
    s) STAGE="${OPTARG}";;
    p) PHASES="${OPTARG}";;
    n) RUN_NUMBER="${OPTARG}";;
    h) help; exit 0;;
    *) usage; exit 1;;
  esac
done

IMPLEMENTED_STAGES="1 2 3"
IMPLEMENTED_PHASES="min eqNVT eqNPT prod NES"

# Check presence of positional arguments
if [ -z "$BASEDIR" ] || [ -z "$TOPDIR" ] || [ -z "$GMX" ] || [ -z "$STAGE" ] || [ -z "$PHASES" ]; then
  usage
  exit 1
fi

# Check option values
mkdir -p "$BASEDIR"
[ -d "$BASEDIR" ] || fail "Directory $BASEDIR does not exist."
[ -d "$TOPDIR" ] || fail "Topology directory $TOPDIR does not exist."
which "$GMX" &> /dev/null || fail "Executable $GMX not found."

VALID_STAGE=0
for stage in $IMPLEMENTED_STAGES; do
  [ "$STAGE" == "$stage" ] && VALID_STAGE=1
done
[ $VALID_STAGE -eq 1 ] || fail "Stage $STAGE is not a valid choice."

for phase in $PHASES; do
  if [[ "$IMPLEMENTED_PHASES" != *"$phase"* ]]; then
    fail "Unknown phase $phase"
  fi
  if [ "$phase" == "min" ]; then
    if [ -z "$CRD" ] || [ ! -e "$CRD" ]; then
      fail "Minimization phase requires a valid coordinate input file."
    fi
  fi
  if [ "$phase" == "NES" ]; then
    if [ -z "$RUN_NUMBER" ] || ! [ "$RUN_NUMBER" -ge 1 ] &> /dev/null; then
      fail "NES phase requires a run number (-n) greater or equal 1"
    fi
    if [ "$STAGE" -ne 2 ] && [ "$STAGE" -ne 3 ]; then
      fail "NES phase is only implemented for stages 2 and 3"
    fi
  fi
done

# Get full paths
BASEDIR=$(readlink -f "$BASEDIR")
TOPDIR=$(readlink -f "$TOPDIR")
# shellcheck disable=SC2236
[ ! -z "$CRD" ] && CRD=$(readlink -f "$CRD")

# Check presence of needed topology files for all stages and phases
for file in system waterRestraintLambdaIndependent waterRestraintLambdaDependent; do
  [ -e "$TOPDIR"/$file.top ] || fail "File $file.top not found in directory $TOPDIR"
done

# Check presence of mdp input files for all stages and phases
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
INPUTDIR=$SCRIPT_PATH/../input
[ -d "$INPUTDIR" ] || fail "Input directory $INPUTDIR not found."
for stage in $IMPLEMENTED_STAGES; do
  [ -e "$INPUTDIR"/stage"$stage".mdp ] || fail "Input file $INPUTDIR/stage$stage.mdp is missing."
done
for phase in $IMPLEMENTED_PHASES; do
  [ -e "$INPUTDIR"/"$phase".mdp ] || fail "Input file $INPUTDIR/$phase.mdp is missing."
done


# Function that runs grompp and mdrun for given directory, mdp file, and gro file
function run_simulation()
{
  local WORKDIR=$1
  local MDP=$2
  local GRO=$3
  local TOP=$4

  STARTDIR=$PWD

  cd "$WORKDIR" || fail "Could not access directory $WORKDIR."

  $GMX grompp -c "$GRO" -p "$TOP" -f "$MDP" || \
    fail "grompp command failed:\n\t$GMX grompp -c \"$GRO\" -p \"$TOP\" -f \"$MDP\"\n\t(wd: $PWD)"
  eval "$GMX mdrun $RUN_PARAMS" || \
    fail "mdrun command failed:\n\t$GMX mdrun \"$RUN_PARAMS\"\n\t(wd: $PWD)"

  cd "$STARTDIR" || fail "Could not access directory $STARTDIR."
}

# Launch simulations for all requested phases
for phase in $PHASES; do
  echo "Running phase $phase ..."

  # Create work directory if it doesn't exist
  if [ "$phase" != "NES" ]; then
    WORKDIR=$BASEDIR/$phase
  else
    WORKDIR=$BASEDIR/$phase/run$RUN_NUMBER
  fi
  mkdir -p "$WORKDIR"

  # Pick input coordinates for phase
  case $phase in
    "min") GRO=$CRD;;
    "eqNVT") GRO=$BASEDIR/min/confout.gro;;
    "eqNPT") GRO=$BASEDIR/eqNVT/confout.gro;;
    "prod") GRO=$BASEDIR/eqNPT/confout.gro;;
    "NES") GRO=$BASEDIR/prod/frames/frame$RUN_NUMBER.gro;;
    *) echo "Unknown phase";;
  esac
  [ -e "$GRO" ] || fail "Phase $phase cannot be run because coordinate file $GRO is missing."

  # Create input file for phase & stage
  MDP=$WORKDIR/input.mdp
  cp "$INPUTDIR"/"$phase".mdp "$MDP" || fail "Error creating input file"
  if [ "$phase" != "NES" ]; then
    cat "$INPUTDIR"/stage"$STAGE".mdp >> "$MDP" || fail "Error creating input file"
  else
    MDP2=$WORKDIR/input2.mdp
    cp "$MDP" "$MDP2"
    cat "$INPUTDIR"/stage"${STAGE}"NES1.mdp >> "$MDP" || fail "Error creating input file"
    cat "$INPUTDIR"/stage"${STAGE}"NES2.mdp >> "$MDP2" || fail "Error creating input file"
  fi

  # Create topology file for phase & stage
  TOP=$WORKDIR/topology.top
  cp "$TOPDIR"/system.top "$TOP" || fail "Error creating topology file"
  if [ "$STAGE" -eq 3 ] || [ "$phase" == "NES" ]; then
    # Stage 3 is running with full restraints, but doesn't need to calculate dH/dL
    # When using NES phase, also stage 2 doesn't need to calculate dH/dL for the restraint
    # If using the lambda dependent restraint, the NES simulations would turn on / off the restraint
    RESTRAINT_FILE=waterRestraintLambdaIndependent.top
  else
    # Stage 1 is running at lambda point that has no restraint
    # Stage 2 is running at lambda point with restraints, calculating dH/dL
    RESTRAINT_FILE=waterRestraintLambdaDependent.top
  fi
  echo >> "$TOP"
  cat "$TOPDIR"/$RESTRAINT_FILE >> "$TOP" || fail "Error creating topology file"

  if [ "$phase" != "NES" ]; then
    run_simulation "$WORKDIR" "$MDP" "$GRO" "$TOP"
  else
    # For NES, we run two simulations, were the second starts from
    # the final configuration of the first
    mkdir -p "$WORKDIR"/1 "$WORKDIR"/2
    run_simulation "$WORKDIR"/1 "$MDP" "$GRO" "$TOP"
    run_simulation "$WORKDIR"/2 "$MDP2" "$WORKDIR"/1/confout.gro "$TOP"
  fi

done
