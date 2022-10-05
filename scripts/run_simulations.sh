#!/bin/bash

IMPLEMENTED_STAGES="1 2 3 3.1 3.2 3.3 3.4 3.5 3.6 3.7 4 5 6 6.1 6.2 6.3 6.4 6.5 6.6 6.7 7"
IMPLEMENTED_PHASES="min eqNVT eqNPT prod NES NES2"

function usage() {
  echo "USAGE: $(basename "$0") [-h] -d dir -t dir [-c crd] -x gmx [-o opts] -s N -p {min,eqNVT,eqNPT,prod,NES,NES2} [-n N] [-m N]"
}

function help() {
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
  # echo -e "\t1: Stage with fully interacting, non-restrained water"
  # echo -e "\t2: Stage with fully interacting, restrained water"
  # echo -e "\t3: Stage with non-interacting, restrained water"
  echo -e "\t$IMPLEMENTED_STAGES"
  echo -e "-p {min,eqNVT,eqNPT,prod,NES,NES2}"
  echo -e "\tThe simulation phase to run, one or more of: $IMPLEMENTED_PHASES"
  echo -e "\tSeparate different phases by spaces, e.g. -p \"min eqNVT eqNPT\""
  echo -e "\tNote: The NES phase is only available for stages 2 through 7"
  echo -e "\tNote: The NES2 phase is only available for stage 5"
  echo -e "-n N\tThe run number for the NES phase"
  echo -e "\tMandatory for the NES phase, ignored otherwise"
  echo -e "-m N\tThe number of steps to run during production, optional"
  echo
}

function fail() {
  echo -e "ERROR: $1"
  exit 1
}

NUM_STEPS=-2
# Get the options
while getopts "d:t:c:x:o:s:p:n:m:h" option; do
  case $option in
    d) BASEDIR="${OPTARG}" ;;
    t) TOPDIR="${OPTARG}" ;;
    c) CRD="${OPTARG}" ;;
    x) GMX="${OPTARG}" ;;
    o) RUN_PARAMS="${OPTARG}" ;;
    s) STAGE="${OPTARG}" ;;
    p) PHASES="${OPTARG}" ;;
    n) RUN_NUMBER="${OPTARG}" ;;
    m) NUM_STEPS="${OPTARG}" ;;
    h)
      help
      exit 0
      ;;
    *)
      usage
      exit 1
      ;;
  esac
done

# Check presence of positional arguments
if [ -z "$BASEDIR" ] || [ -z "$TOPDIR" ] || [ -z "$GMX" ] || [ -z "$STAGE" ] || [ -z "$PHASES" ]; then
  usage
  exit 1
fi

# Check option values
mkdir -p "$BASEDIR"
[ -d "$BASEDIR" ] || fail "Directory $BASEDIR does not exist."
[ -d "$TOPDIR" ] || fail "Topology directory $TOPDIR does not exist."
which "$GMX" &>/dev/null || fail "Executable $GMX not found."

VALID_STAGE=0
for stage in $IMPLEMENTED_STAGES; do
  [ "$STAGE" == "$stage" ] && VALID_STAGE=1
done
[ $VALID_STAGE -eq 1 ] || fail "Stage $STAGE is not a valid choice."
LAMBDA_WINDOW=${STAGE:2:1}
[ -n "$LAMBDA_WINDOW" ] && STAGE=${STAGE:0:1}

for phase in $PHASES; do
  if [[ "$IMPLEMENTED_PHASES" != *"$phase"* ]]; then
    fail "Unknown phase $phase"
  fi
  if [ "$phase" == "min" ]; then
    if [ -z "$CRD" ] || [ ! -e "$CRD" ]; then
      fail "Minimization phase requires a valid coordinate input file."
    fi
  fi
  if [ "$phase" == "NES" ] || [ "$phase" == "NES2" ]; then
    if [ -z "$RUN_NUMBER" ] || ! [ "$RUN_NUMBER" -ge 1 ] &>/dev/null; then
      fail "NES phase requires a run number (-n) greater or equal 1"
    fi
  fi
  if [ "$phase" == "NES2" ]; then
    [ "$STAGE" == "5" ] || fail "NES2 phase is only implemented for stage 5"
  fi
done

# Get full paths
BASEDIR=$(readlink -f "$BASEDIR")
TOPDIR=$(readlink -f "$TOPDIR")
# shellcheck disable=SC2236
[ ! -z "$CRD" ] && CRD=$(readlink -f "$CRD")

# Check presence of needed topology files for all stages and phases
[ -e "$TOPDIR"/system.top ] || fail "File system.top not found in directory $TOPDIR"

# Check presence of mdp input files for all stages and phases
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
INPUTDIR=$SCRIPT_PATH/../input
[ -d "$INPUTDIR" ] || fail "Input directory $INPUTDIR not found."
for stage in $IMPLEMENTED_STAGES; do
  [ -n "${stage:2:1}" ] && continue
  [ -e "$INPUTDIR"/stage"$stage".mdp ] || fail "Input file $INPUTDIR/stage$stage.mdp is missing."
done
for phase in $IMPLEMENTED_PHASES; do
  [ -e "$INPUTDIR"/"$phase".mdp ] || fail "Input file $INPUTDIR/$phase.mdp is missing."
done

# Function that takes care of multiple `define = ...` statements in mdp files
function unify_define_statement() {
  local MDP=$1
  # Check whether there's more than one define statement
  num_define_statements=$(egrep '^define *=' $MDP | wc -l)
  [ "$num_define_statements" -lt 2 ] && return
  # Collect existing define statements
  define_statements=""
  while read -r line; do
    define_statements="${define_statements} ${line##*=}"
  done < <(grep -E '^define *=' "$MDP")
  # Remove existing define statements
  sed -i.bak '/^define *=/d' "$MDP"
  rm "$MDP.bak"
  # Write consolidated define statement
  echo ""
  echo "; Consolidated define statement" >> "$MDP"
  echo "define = $define_statements" >> "$MDP"
}

# Function that runs grompp and mdrun for given directory, mdp file, and gro file
function run_simulation() {
  local WORKDIR=$1
  local MDP=$2
  local GRO=$3
  local TOP=$4
  local WARNINGS=$5
  local POSRES=$6
  local NSTEPS=$7

  STARTDIR=$PWD

  cd "$WORKDIR" || fail "Could not access directory $WORKDIR."

  if [ -n "$POSRES" ]; then
    [ -e "$POSRES" ] || fail "Positional restraint file $POSRES is missing."
    POSRES="-r $POSRES"
  fi

  GROMPP_CMD="$GMX grompp -c $GRO -p $TOP -f $MDP -maxwarn $WARNINGS $POSRES"
  eval "$GROMPP_CMD" || fail "grompp command failed:\n\t$GROMPP_CMD\n\t(wd: $PWD)"
  MDRUN_CMD="$GMX mdrun -nsteps $NSTEPS $RUN_PARAMS"
  eval "$MDRUN_CMD" ||
    fail "mdrun command failed:\n\t$MDRUN_CMD\n\t(wd: $PWD)"

  cd "$STARTDIR" || fail "Could not access directory $STARTDIR."
}

# Launch simulations for all requested phases
for phase in $PHASES; do
  echo "Running phase $phase ..."

  # Create work directory if it doesn't exist
  if [ "$phase" != "NES" ] && [ "$phase" != "NES2" ]; then
    WORKDIR=$BASEDIR/$phase
  else
    WORKDIR=$BASEDIR/$phase/run$RUN_NUMBER
  fi
  mkdir -p "$WORKDIR"

  # Pick input coordinates for phase
  case $phase in
    "min") GRO=$CRD ;;
    "eqNVT") GRO=$BASEDIR/min/confout.gro ;;
    "eqNPT") GRO=$BASEDIR/eqNVT/confout.gro ;;
    "prod") GRO=$BASEDIR/eqNPT/confout.gro ;;
    "NES") GRO=$BASEDIR/prod/frames/frame$RUN_NUMBER.gro ;;
    "NES2") GRO=$BASEDIR/prod/frames/frame$RUN_NUMBER.gro ;;
    *) echo "Unknown phase" ;;
  esac
  [ -e "$GRO" ] || fail "Phase $phase cannot be run because coordinate file $GRO is missing."

  # Create input file for phase & stage
  MDP=$WORKDIR/input.mdp
  cp "$INPUTDIR"/"$phase".mdp "$MDP" || fail "Error creating input file"
  if [ "$phase" != "NES" ] && [ "$phase" != "NES2" ]; then
    cat "$INPUTDIR"/stage"$STAGE".mdp >>"$MDP" || fail "Error creating input file"
  elif [ "$phase" == "NES" ]; then
    cat "$INPUTDIR"/stage"${STAGE}"NES.mdp >>"$MDP" || fail "Error creating input file"
  elif [ "$phase" == "NES2" ]; then
    cat "$INPUTDIR"/stage"${STAGE}"NES2.mdp >>"$MDP" || fail "Error creating input file"
  fi
  if [ -n "$LAMBDA_WINDOW" ]; then
    perl -pi -e "s/init_lambda_state        = 2/init_lambda_state        = $((LAMBDA_WINDOW+2))/" $MDP
  fi
  if [ -e "$TOPDIR"/system.mdp ] && [ "$phase" != "NES2" ]; then
    cat "$TOPDIR"/system.mdp >>"$MDP" || fail "Error creating input file"
  fi
  if [ "$phase" == "NES2" ]; then
    cat "$TOPDIR"/systemNES2.mdp >>"$MDP" || fail "Error creating input file"
  fi

  unify_define_statement "$MDP"

  # Create topology file for phase & stage
  TOP="$TOPDIR"/system.top

  [ "$phase" = "min" ] && WARNINGS=2 || WARNINGS=0

  POSRES="$TOPDIR/minimized.gro"
  if [ "$phase" = "eqNVT" ] || [ "$phase" = "eqNPT" ]; then
    POSRES="$BASEDIR/min/confout.gro"
  fi

  [ "$phase" = "prod" ] && NSTEPS=$NUM_STEPS || NSTEPS=-2

  run_simulation "$WORKDIR" "$MDP" "$GRO" "$TOP" "$WARNINGS" "$POSRES" "$NSTEPS"

done
