# waterNES
This repository contains workflows to calculate relative free energies
using non-equilibrium switching for buried water molecules. These
workflows operate on thermodynamic cycles like these:

![Thermodynamic cycle](https://github.com/MobleyLab/waterNES/blob/main/docs/ThermodynamicCycle.pdf?raw=true)

The purple cycle includes additional restraints for the protein, and
is preferred when it is likely that the binding site might collapse
once the water is removed. In all other cases, the simpler blue cycle
can be used.

### Cycle branches
* **A**, the target branch: We will find its value by summing all
  other branches.
* **B**, water restraints correction for interacting water: Calculated
  from equilibrium simulations of its end points (stages 1 and 2,
  simulations of the protein interacting with the unrestrained and
  restrained water, respectively).
* **C**, water coupling / decoupling: Calculated using multiple NES
  simulations from each end point to the other (stages 2 and 3,
  simulations of the protein with and without interactions with the
  restrained water). The starting structures are obtained from long
  equilibrium simulations at the end points.
* **D**, water restraints correction for non-interacting water: This
  value can be calculated analytically, no simulations required.
* **E**, transfer of non-interacting water to bulk water: This branch
  does not involve any free energy change, $\Delta G_{E} = 0$.
* **F**, turn on interactions of water molecule in bulk water: This
  process is only depending on the water model, so it can be simulated
  once and reused for every system which uses the same water model.
  This is currently not part of the workflow.
  
### Stages
* **Stage 1**: This simulates the protein complex with an
  unrestained, fully interacting water molecule.  Simulations include
  minimization, equilibration in NVT and NPT, production simulation in
  NPT. During the production run, the free energy with respect to the
  water restraint (which is turned off) is calculated.
* **Stage 2**: This simulates the protein complex with a restrained,
  fully interacting water molecule. Simulations include minimization,
  equilibration in NVT and NPT, production simulation in NPT. During
  the production run, the free energy with respect to the water
  restraint (which is turned on) is calculated. Characteristic
  structures from the production simulation are used to start NES
  simulations (transition stage 2 to 3).
* **Stage 3**: This simulates the protein complex with a restrained,
  non-interacting water molecule. Simulations include minimization,
  equilibration in NVT and NPT, production simulation in
  NPT. Characteristic structures from the production simulation are
  used to start NES simulations (transition stage 3 to 2).
* **Stages 4 and 5**: Coming soon.

### Systems
Currently, this repo contains input files for five proteins, each in
combination with one to three different ligands. The protein-ligand
complexes are named after their PDB ID. They include

* HIV-1 protease: 1HPX, 1EC0, 1EBW
* Trypsin: 1AZ8, 1C5T, 1GI1
* Scytalone dehydratase[1]: 3STD, 4STD, 7STD
* Factor Xa: 1EZQ, 1LPG, 1F0S
* BPTI: 5PTI

[1]: 3STD has a different water site than 4STD and 7STD

## Usage
### Run simulations
All simulations can be run using the `bash` file
`scripts/run_simulations.sh`. There is an integrated documentation
which can be called by

```
$ bash scripts/run_simulations.sh -h
Run Water NES simulations.

USAGE: run_water_nes [-h] -d dir -t dir [-c crd] -x gmx [-o opts] -s N -p {min,eqNVT,eqNPT,prod}
Options:
-h	Print this help and exit
-d dir	The base working directory for the system
-t top	The topology directory including topology files
-c crd	The input coordinates, for minimization only
-x gmx	The GROMACS executable to use
-o opts	The run options, e.g. "-nt 16 -ntomp 4", optional
-s N	The stage to run, one of:
	1: Stage with fully interacting, non-restrained water
	2: Stage with fully interacting, restrained water
	3: Stage with non-interacting, restrained water
-p {min,eqNVT,eqNPT,prod}
	The simulation phase to run, one or more of: min, eqNVT, eqNPT, prod
	Separate different phases by spaces, e.g. -p "min eqNVT eqNPT".
```

In the above `USAGE` string, flags in brackets (`[]`) are optional,
options in braces (`{}`) denote possible values, and all other flags
are mandatory. Note that the coordinates (`-c` flag) is mandatory
for the minimization phase, and ignored for all other phases.

The directory (`-d`) is the base directory. All simulations input and
output is going to be written in subdirectories of this base
directory. The only requirement for this directory is full write
access. All files within the directory might be overwritten.

The topology directory (`-t`) needs to contain topology files for the
system to be simulated. Required files are the main topology named as
`system.top`, as well as files containing $\lambda$-dependent and
$\lambda$-independent position restraints for the water, named
`waterConstraintLambdaIndependent.top` and
`waterConstraintLambdaDependent.top`, respectively.

The coordinate file (`-c`) needs to contain starting coordinates for
the system matching with the topology found in the topology directory.
This is only needed for the mimization phase.

The GROMACS executable to use (`-x`) can either be a path to an
executable, or simply `gmx` if it is in the path.

The run options (`-o`) are only used for the call to `gmx mdrun`, and
allow to set parallelization details or other run time options
accepted by `gmx mdrun`.

The stage (`-s`) denotes one of the stages of the thermodynamic cycle
as defined above.

The simulation phase (`-p`) defines the simulation to be run, one of
`min` (minimization of the input configuration), `eqNVT`
(thermalization at constant volume), `eqNPT` (equilibration at
constant temperature and pressure), and `prod` (production simulation
in NPT). Each phase builds on the previous phase. Several phases can
be combined.

### Analysis
Coming soon.

### Input files
#### GROMACS mdp files
The input parameters for each of the phases (min / eqNVT / eqNPT /
prod) are almost identical between the different stages of the
thermodynamics cycle. The common parameters are stored in
`input/{min,eqNVT,eqNPT,prod}.mdp`. Any changes to the simulations of
the different phases can be done in these files.

The difference between the stages are determined by the free energy
settings of the parameter file. These are stored in
`input/stage{1,2,3}.mdp`. Any changes to the free energy settings
defining the different phases can be done in these files.

The `run_simulations.sh` script combines the phase-defining and the
stage-defining mdp files to create the run file used by `gmx grompp`.

#### GROMACS topology files
The topology files for the different systems mentioned are stored in
folders named by their PDB ID in `systems/`. The base topologies
(named `system.top`) do not contain the intramolecular interaction
which restrains the water to the protein, because these restraints
vary based on the stage in the thermodynamic cycle. These interactions
are stored in separate files in the same directories named
`waterConstraintLambdaIndependent.top` and
`waterConstraintLambdaDependent.top`, respectively.

#### Configuration files
The input configurations (named `systems.gro`) used for the
minimization are stored in the same folders as their respective
topology files.

## Tests
### `tests/test_run_simulations.sh`
This script tests the `scripts/run_simulations.sh` file by making sure
that the typical usage (including calls to GROMACS) runs without
errors. The test is invoked by calling `bash
tests/test_run_simulations.sh` from the root of the repository. The
test script expects `gmx` to be in the path.

Github actions runs this test on every push and every PR to `main`.
