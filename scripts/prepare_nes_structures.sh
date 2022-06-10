#!/bin/bash

function usage()
{
  echo "USAGE: $(basename $0) [-h] -d dir -x gmx -n num"
}

function help()
{
  echo "Prepare structures for NES simulations"
  echo
  usage
  echo "Options:"
  echo -e "-h\tPrint this help and exit"
  echo -e "-d dir\tThe base working directory of the stage"
  echo -e "\tThis directory is expected to have a subdirectory prod/"
  echo -e "\twith the results of the production simulation of that stage"
  echo -e "-x gmx\tThe GROMACS executable to use"
  echo -e "-n num\tThe number of structures to prepare"
  echo
}

function fail()
{
  echo -e "ERROR: $1"
  exit 1
}

# Get the options
while getopts "hd:x:n:" option; do
  case $option in
    d) BASEDIR="${OPTARG}";;
    x) GMX="${OPTARG}";;
    n) NUM_STRUCTURES="${OPTARG}";;
    h) help; exit 0;;
    *) usage; exit 1;;
  esac
done

# Check presence of positional arguments
if [ -z "$BASEDIR" ] || [ -z "$GMX" ] || [ -z "$NUM_STRUCTURES" ]; then
  usage
  exit 1
fi

# Check option values
[ -d "$BASEDIR" ] || fail "Directory $BASEDIR does not exist."
which $GMX &> /dev/null || fail "Executable $GMX not found."
[ $NUM_STRUCTURES -gt 0 ] &> /dev/null || fail "Number of structures $NUM_STRUCTURES is not a valid number."

# Make sure all relevant files and folders are present
WORKDIR=$BASEDIR/prod
[ -d $WORKDIR ] || fail "Production simulation folder $WORKDIR not found."
TRJ=$WORKDIR/traj.trr
[ -e $TRJ ] || fail "Production trajectory $TRJ not found."
MDP=$WORKDIR/mdout.mdp
[ -e $MDP ] || fail "Production mdp file $MDP not found."

# Find the frequency at which we want to extract structures, do some sanity checks
WRITE_FREQUENCY=$(grep "nstxout " $MDP | awk '{print $NF;}')
NUM_STEPS=$(grep "nsteps " $MDP | awk '{print $NF;}')
NUM_FRAMES=$(($NUM_STEPS/$WRITE_FREQUENCY))
echo "INFO: Production simulation with nsteps = $NUM_STEPS, nstxout = $WRITE_FREQUENCY, generated $NUM_FRAMES frames."

[ $NUM_FRAMES -ge $NUM_STRUCTURES ] || \
  fail "Production simulation generated only $NUM_FRAMES frame, cannot create $NUM_STRUCTURES structures."
[ $(($NUM_FRAMES%$NUM_STRUCTURES)) -eq 0 ] || \
  echo "WARNING: Production simulation generated $NUM_FRAMES, which is not a multiple of the number of structures ($NUM_STRUCTURES)."

# Go to simulation directory
STARTDIR=$PWD
cd $WORKDIR || fail "Could not access directory $WORKDIR."
mkdir -p frames

# Extract frames
echo "System" | $GMX trjconv -f traj.trr -o frames/frame.gro -sep -skip $(($NUM_FRAMES/$NUM_STRUCTURES)) -ur compact -pbc mol || \
  fail "trjconv command failed:\n\t$GMX trjconv -f $TRJ -o frame.gro -sep -skip $(($NUM_FRAMES/$NUM_STRUCTURES)) -ur compact -pbc mol"
# Ignore the zero frame
rm -f frames/frame0.gro

# Return to initial directory & exit
cd $STARTDIR || fail "Could not access directory $STARTDIR."
FRAME_NUMBERS=$([ $NUM_STRUCTURES -eq 1 ] && echo 1 || echo "{1-$NUM_STRUCTURES}")
echo "Successfully created $NUM_STRUCTURES frames, $WORKDIR/frames/frame$FRAME_NUMBERS.gro"
