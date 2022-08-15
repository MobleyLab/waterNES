#!/bin/bash

IMPLEMENTED_STAGES="1 1a 1b 1c 1d 2 2a 3 3a 4a 4b 5a 5b 6 6a 7 7a 8b 8c 8d"
IMPLEMENTED_PHASES="min eqNVT eqNPT prod NES"

function usage() {
  echo "USAGE: $(basename "$0") [-h] -d dir -t dir [-c crd] -x gmx [-o opts] -s N -p {min,eqNVT,eqNPT,prod,NES} [-n N] [-m N]"
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
  echo -e "-p {min,eqNVT,eqNPT,prod}" #,NES}"
  echo -e "\tThe simulation phase to run, one or more of: $IMPLEMENTED_PHASES"
  echo -e "\tSeparate different phases by spaces, e.g. -p \"min eqNVT eqNPT\""
  # echo -e "\tNote: The NES phase is only available for stages 2 and 3"
  # echo -e "-n N\tThe run number for the NES phase"
  # echo -e "\tMandatory for the NES phase, ignored otherwise"
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
    if [ -z "$RUN_NUMBER" ] || ! [ "$RUN_NUMBER" -ge 1 ] &>/dev/null; then
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
[ -e "$TOPDIR"/system.top ] || fail "File system.top not found in directory $TOPDIR"

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
  eval "$GROMPP_CMD" ||
    fail "grompp command failed:\n\t$GROMPP_CMD\n\t(wd: $PWD)"
  MDRUN_CMD="$GMX mdrun -nsteps $NSTEPS $RUN_PARAMS"
  eval "$MDRUN_CMD" ||
    fail "mdrun command failed:\n\t$MDRUN_CMD\n\t(wd: $PWD)"

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
    "min") GRO=$CRD ;;
    "eqNVT") GRO=$BASEDIR/min/confout.gro ;;
    "eqNPT") GRO=$BASEDIR/eqNVT/confout.gro ;;
    "prod") GRO=$BASEDIR/eqNPT/confout.gro ;;
    "NES") GRO=$BASEDIR/prod/frames/frame$RUN_NUMBER.gro ;;
    *) echo "Unknown phase" ;;
  esac
  [ -e "$GRO" ] || fail "Phase $phase cannot be run because coordinate file $GRO is missing."

  # Create input file for phase & stage
  MDP=$WORKDIR/input.mdp
  cp "$INPUTDIR"/"$phase".mdp "$MDP" || fail "Error creating input file"
  if [ "$phase" != "NES" ]; then
    cat "$INPUTDIR"/stage"$STAGE".mdp >>"$MDP" || fail "Error creating input file"
  else
    MDP2=$WORKDIR/input2.mdp
    cp "$MDP" "$MDP2"
    cat "$INPUTDIR"/stage"${STAGE}"NES1.mdp >>"$MDP" || fail "Error creating input file"
    cat "$INPUTDIR"/stage"${STAGE}"NES2.mdp >>"$MDP2" || fail "Error creating input file"
    unify_define_statement "$MDP2"
  fi

  unify_define_statement "$MDP"

  stage_modifier=${STAGE:1:1}
  if [ "$phase" = "prod" ]; then
    if [ -n "$stage_modifier" ]; then
      cat "$INPUTDIR"/stageNxprod.mdp >>"$MDP" || fail "Error creating input file"
    else
      cat "$INPUTDIR"/stageNprod.mdp >>"$MDP" || fail "Error creating input file"
    fi
  fi

  # Create topology file for phase & stage
  TOP="$TOPDIR"/system.top

  [ "$phase" = "min" ] && WARNINGS=1 || WARNINGS=0

  POSRES=""
  if [ -n "$stage_modifier" ]; then
    if [ "$phase" != "min" ]; then
      POSRES="$TOPDIR/min/restraint.gro"
    fi
  fi
  if [ "$phase" = "eqNVT" ] || [ "$phase" = "eqNPT" ]; then
    POSRES="$BASEDIR/min/confout.gro"
  fi

  [ "$phase" = "prod" ] && NSTEPS=$NUM_STEPS || NSTEPS=-2

  if [ "$phase" != "NES" ]; then
    run_simulation "$WORKDIR" "$MDP" "$GRO" "$TOP" "$WARNINGS" "$POSRES" "$NSTEPS"
  else
    # For NES, we run two simulations, where the second starts from
    # the final configuration of the first
    mkdir -p "$WORKDIR"/1 "$WORKDIR"/2
    run_simulation "$WORKDIR"/1 "$MDP" "$GRO" "$TOP" "$WARNINGS" "$POSRES" "$NSTEPS"
    run_simulation "$WORKDIR"/2 "$MDP2" "$WORKDIR"/1/confout.gro "$TOP" "$WARNINGS" "$POSRES" "$NSTEPS"
  fi

done
