# waterNES

**A reproducible, modular framework for free energy calculations in protein-ligand systems involving trapped/buried water molecules using non-equilibrium switching (NES) and equilibrium free energy simulations.**

The methods implemented here can perform the following free energy calculations:
- Relative binding free energy (RBFE) calculation between two protein-ligand complexes involving trapped/buried water molecules in ligand binding.
  - *This has two versions: a fully NES version (Check thermodynamic cycle 1 below), and an NES + equilibrium free energy calculation version (Check thermodynamic cycle 2 below.).*
- Absolute binding free energy (ABFE) calculation of trapped/buried water molecules in proteins or protein-ligand complexes.
  - *This is performed using NES + equilibrium free energy calculation only (Check thermodynamic cycle 3 below).*
- Fullerene free energy calculations. 
  - *NES free energy calculation of displacing a water molecule from a C<sub>90</sub> fullerene cavity, performed using  free energy calculations.*
  - *Free energy calculation of solvent water equilibration inside a C<sub>90</sub> fullerene cavity using Hamiltonia replica exchange (HREX) simulation.*

### Associated publications:

1. *[Fullerene free energy calculations]* A Quadrupolar Fullerene Model System for Benchmarking Enhanced Sampling of Trapped Waters in Free Energy Calculations; Swapnil Wagle and David L. Mobley; *J. Phys. Chem. B* **2026**, 130, 2869-2882 ([Fullerene Paper](https://pubs.acs.org/doi/full/10.1021/acs.jpcb.5c08189))

2. *[RBFE trapped/buried waters- fully NES version]* Advancing Binding Affinity Calculations: A Non-Equilibrium Simulations Approach for Calculation of Relative Binding Free Energies in Systems with Trapped Waters; Swapnil Wagle, Christopher I. Bayly and L. David Mobley; *J. Chem. Theory Comput.* **2025**, 21, 7593-7604 ([RBFE NES Paper](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00758))

3. *[RBFE trapped/buried waters- NES + equilibrium free energy version,  ABFE of trapped/buried waters]* Leveraging a Separation of States Method for Relative Binding Free Energy Calculations in Systems with Trapped Waters;
Swapnil Wagle, Pascal T. Merz, Yunhui Ge, Christopher I. Bayly and David L. Mobley; *J. Chem. Theory Comput.* **2024**, 20, 11013−11031 ([ABFE-RBFE Paper](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01145))


## :rocket: Quick start

### RBFE calculation using NES simulations



The following image shows thermodynamic cycle 1 to calculate RBFE between a ligand pair involving trapped water molecules in ligand binding, implementing a non-equilibrium switching (NES)-based workflow.
<img src="https://github.com/MobleyLab/waterNES/blob/main/docs/NES-Total-Chemdraw.png?raw=true" style="width:70%;" alt="Thermodynamic cycle 1">

Here, ligand A is shown in green and ligand B is shown in purple. Both ligands are bound to the protein (shown in yellow). A trapped water is shown in red, while a decoupled trapped water is shown in pale red.
A black cross on the trapped water represents harmonic restraint, a dashed circle around the trapped water represents solvent repulsion potential applied at the binding site of the trapped water.
NES-Total, H and J consitute the complex leg of the RBFE calculation. NES-Solvent represents the solvent leg of the RBFE calculation. 

In the thermodynamic cycle, we simulate the complex and the solvent legs of the RBFE calculation. For the complex leg, we will simulate Stage 1 and Stage 4 for 6 ns. From each of the two trajectories, we extract 100 frames and will launch the NES switches. Each of the NES switches will run for 70 ps (10 ps, 50 ps, 10 ps for NES1, NES2 and NES3, respectively; one NES1, NES2 and NES3 switch constitute one NES-Total switch). For the solvent leg, we will simulate each of the two ligands in solvent for 6 ns and will extract 100 frames from each trajectory. From each frame, we will launch the NES-Solvent switch for 50 ps.

Finally, we will calculate the RBFE between the two ligands from these NES trajectories, using the ```analysis``` module. 

RBFE = $ΔG_{NES-Total} + \Delta G_H + \Delta G_J -  \Delta G_{NES-Solvent}$

## Usage

Here, we are showing how the workflow is used by giving an example of RBFE calculation between ligands C3d and C5d (check Figure 3 of the [RBFE NES Paper](https://pubs.acs.org/doi/full/10.1021/acs.jctc.5c00758) for their structures).   

To start, run the following commands to download the repository and create a conda environment with the necessary packages.

```
git clone git@github.com:MobleyLab/waterNES.git
cd waterNES
conda env create -f environment.yml
conda activate waterNES
```
To run the RBFE calculation, we need two files: a structure of a dual topology protein-ligand complex (```minimized.gro```) and its topology file (```system.top```). Both these files are given for the example case in the folder ```examples/rbfe_nes/sd_c3d_c5d/complex_leg```. Similarly, to run the solvent leg of the calculation, the structure of topology files are saved in ```examples/rbfe_nes/sd_c3d_c5d/solvent_leg```.

Now, we will run the ```run.py``` file that will start the end state simulations of Stage 1 and Stage 4, and the end state simulations of the solvent leg. In each end state run, first the system will be energy minimized, and then will be equilibrated in an NVT and then in an NPT ensemble. If ```posre.itp``` and/or ```posre_ligand.itp``` files are present in the directory (directory structure shown below), then they will be used in the NVT and NPT equilibrations. Finally, 6 ns of production run will be produced starting from the NPT equilibrated structure. The 6 ns trajectory will be divided into 100 structures, each of which will be used to launch NES runs (NES1 -> NES2 -> NES3 for the complex leg, and NES-Solvent for the solvent run). 

```
python run.py rbfe examples/rbfe_nes/config.yaml
```
The run.py will read the ```minimized.gro``` and the ```system.top``` files and will generate the output according to the following directoty structure.

For this example run, we have already performed the simulations, so you don't have to perform the whole runs again (takes around 1.5 hour when run on 20 CPUs and 1 GPU).
The important files are the ```dhdl.xvg``` files generated by the NES switches. We will use these files into the analyze our simulations and calculate the RBFE.

To calculate the RBFE for this example run, we will run the following command:

```
python analysis/analyze.py rbfe examples/rbfe_nes/config.yaml
```
The ```analyze.py``` script takes in all the ```dhdl.xvg``` as input and calculates the free energies $ΔG_{NES-Total}$ and $ΔG_{NES-Solvent}$ using MBAR calculations.

The expected output of the command is:

```RBFE = -1.98 +- 0.1 kcal/mol```

The following graphs will also be generated, and will be saved in the ```examples/rbfe_nes/sd_c3d_c5d/complex_leg``` and ```examples/rbfe_nes/sd_c3d_c5d/solvent_leg``` folders. These graphs show the histograms of the work values in the forward and reverse NES switches. For 

Thermodynamic cycle 2 to calculate RBFEs between ligands that bind to target protein with different numbers of
trapped water molecules. The RBFE is calculated using the following thermodynamic cycle
![Thermodynamic cycle](https://github.com/MobleyLab/waterNES/blob/main/docs/rbfe_cycle_short.png?raw=true)

In addition to the description for the thermodynamic cycle 1, te black cross on the protein represents position restraints on the binding site. In this figure, only the complex leg of the RBFE calculation is shown.
The solvent leg of the RBFE calculation is same as shown for the thermodynamic cycle for 1. 

Thermodynamic cycle 3 to calculate absolute binding free energies (ABFEs) of trapped waters in proteins/protein-ligand complexes using
the following thermodynamic cycle
![Thermodynamic cycle](https://github.com/MobleyLab/waterNES/blob/main/docs/abfe_cycle_full.png?raw=true)

Here, the absolute binding free energy of a trapped water is calculated, i.e., the free energy to displace a water from a protein's binding site. 

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
  unrestrained, fully interacting water molecule.  Simulations include
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
which can be called by running the script with the `-h` flag:

```shell
$ bash scripts/run_simulations.sh -h
Run Water NES simulations.

USAGE: run_simulations.sh [-h] -d dir -t dir [-c crd] -x gmx [-o opts] -s N -p {min,eqNVT,eqNPT,prod,NES} [-n N]
Options:
-h      Print this help and exit
-d dir  The base working directory for the stage
-t top  The topology directory including topology files
-c crd  The input coordinates
        Mandatory for the minimization phase, ignored otherwise
-x gmx  The GROMACS executable to use
-o opts The run options, e.g. "-nt 16 -ntomp 4", optional
-s N    The stage to run, one of:
        1: Stage with fully interacting, non-restrained water
        2: Stage with fully interacting, restrained water
        3: Stage with non-interacting, restrained water
-p {min,eqNVT,eqNPT,prod,NES}
        The simulation phase to run, one or more of: min, eqNVT, eqNPT, prod, NES
        Separate different phases by spaces, e.g. -p "min eqNVT eqNPT"
        Note: The NES phase is only available for stages 2 and 3
-n N    The run number for the NES phase
        Mandatory for the NES phase, ignored otherwise
```

In the above `USAGE` string, flags in brackets (`[]`) are optional,
options in braces (`{}`) denote possible values, and all other flags
are mandatory.

The directory (`-d`) is the base directory. All simulations input and
output is going to be written in subdirectories of this base
directory. The only requirement for this directory is full write
access. All files within the directory might be overwritten.

The topology directory (`-t`) needs to contain topology files for the
system to be simulated. Required files are the main topology named as
`system.top`, as well as files containing $\lambda$-dependent and
$\lambda$-independent position restraints for the water, named
`waterRestraintLambdaIndependent.top` and
`waterRestraintLambdaDependent.top`, respectively.

The coordinate file (`-c`) needs to contain starting coordinates for
the system matching with the topology found in the topology directory.
This is only needed for the minimization phase.

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
constant temperature and pressure), `prod` (production simulation in
NPT), and `NES` (non-equilibrium transition simulations, only for
stages 2 and 3). Each phase builds on the previous phase. Several
phases can be combined. Note that the `NES` phase requires that the
starting structures were prepared from the result of the `prod` phase,
so it cannot be directly chained.

For the NES phase, the run number needs to be given with the flag
`-n`. Each run starts from a different configuration, extracted from
the production simulation with the `scripts/prepare_nes_structures.sh`
script. The `-n` flag is ignored for all other phases.

### Prepare configurations for NES runs
NES simulations need to be run from multiple starting structures
representative of the ensemble of configurations in the end states. To
achieve this, we extract configurations from the long NPT production
simulations in stages 2 and 3. This operation can be performed using
the `scripts/prepare_nes_structures.sh` script. Again, a help function
describes the functionality:

```shell
Prepare structures for NES simulations

USAGE: prepare_nes_structures.sh [-h] -d dir -x gmx -n num
Options:
-h      Print this help and exit
-d dir  The base working directory of the stage
        This directory is expected to have a subdirectory prod/
        with the results of the production simulation of that stage
-x gmx  The GROMACS executable to use
-n num  The number of structures to prepare
```

In the above `USAGE` string, flags in brackets (`[]`) are optional,
while all other flags are mandatory.

The directory (`-d`) is the base directory of the stage in which all
simulation input and output was written, and should be the same
directory that was used for the `-d` flag of the `run_simulations.sh`
script. Specifically, the script is looking for the mdp file and the
trajectory file in the `prod/` folder under this base directory.

The GROMACS executable to use (`-x`) can either be a path to an
executable, or simply `gmx` if it is in the path.

Finally, `-n` denotes the number of structures to prepare. This
defines how many structures are extracted from the production
simulation for further use.

### Analysis
#### Edge B
Analysis of the equilibrium endpoint free energy difference is implemented
in a Python function. See the following simple example for its usage:
```python
import pathlib

from water_nes.analysis.endpoint_free_energy import calculate_endpoint_free_energy

# Define input directory - here we use the test input files
input_directory = pathlib.PurePath('tests/endpoint_input_files')

# The GROMACS free energy xvg files for the end states (lambda = 0 and
# lambda = 1) are stored in files `lambda0.xvg` and `lambda1.xvg`.
# In the following example, we calculate the endpoint free energy in
# kcal/mol, ignoring the initial 1ns and any part of the trajectory
# after 15ns.
free_energy_estimate = calculate_endpoint_free_energy(
    file_lambda_0=input_directory.joinpath("lambda0.xvg"),
    file_lambda_1=input_directory.joinpath("lambda1.xvg"),
    start_time=1000,
    end_time=15000,
    output_units="kcal/mol",
)

# Print free energy estimate
print(f"Free energy estimate: "
      f"{free_energy_estimate.value:.2f} += {free_energy_estimate.error:.2f} "
      f"{free_energy_estimate.units}")
```

The last line prints the estimate obtained by the `calculate_endpoint_free_energy`
function. For the test files, it will look like this:
```shell
Free energy estimate: 2.10 += 0.28 kcal/mol
```

#### Edge C
Analysis of the non-equilibrium switching calculations is implemented in
a Python function. Here's a minimal example:
```python
import pathlib

from water_nes.analysis.nes_free_energy import calculate_nes_free_energy

# Define input directory - here we use the test input files
input_directory = pathlib.PurePath('tests/nes_input_files')

# Create lists of input files. The transition between stages A and B happens in
# two steps, by first turning off the Coulomb interactions, then the vdW
# interactions, or vice-versa.
xvg_forward_coulomb = [
    input_directory.joinpath(f"transition_A2B_coul_{n}.xvg") for n in range(1, 11)
]
xvg_forward_vdw = [
    input_directory.joinpath(f"transition_A2B_vdw_{n}.xvg") for n in range(1, 11)
]
xvg_backward_coulomb = [
    input_directory.joinpath(f"transition_B2A_coul_{n}.xvg") for n in range(1, 11)
]
xvg_backward_vdw = [
    input_directory.joinpath(f"transition_B2A_vdw_{n}.xvg") for n in range(1, 11)
]

# Get free energy estimate in kcal/mol
free_energy_estimate = calculate_nes_free_energy(
    xvg_files_forward_transition=[xvg_forward_coulomb, xvg_forward_vdw],
    xvg_files_backward_transition=[xvg_backward_vdw, xvg_backward_coulomb],
    temperature=298.15,
    output_units="kcal/mol",
    bootstrapping_repeats=0,
)

# Print free energy estimate
print(f"Free energy estimate: "
      f"{free_energy_estimate.value:.2f} += {free_energy_estimate.error:.2f} "
      f"{free_energy_estimate.units}")
```

The last line prints the estimate obtained by the `calculate_nes_free_energy`
function. For the test files, it will look like this:

```shell
Free energy estimate: 10.41 += 0.62 kcal/mol
```

This is the free energy estimate for edge C, the transition from stage 2 to stage 3.

### Input files
#### GROMACS mdp files
The input parameters for each of the phases (min / eqNVT / eqNPT /
prod / NES) are almost identical between the different stages of the
thermodynamics cycle. The common parameters are stored in
`input/{min,eqNVT,eqNPT,prod,nes}.mdp`. Any changes to the simulations
of the different phases can be done in these files.

The difference between the stages are determined by the free energy
settings of the parameter file. These are stored in
`input/stage{1,2,3}.mdp`. Any changes to the free energy settings
defining the different phases can be done in these files.

For the NES phase, the free energy settings are different from the
other stages, and there are two independent simulations (coupling or
decoupling the vdW and the Coulomb interactions separately). As a
result, there are two additional parameter files for stages 2 and 3,
`input/stage{2,3}NES{1,2}.mdp`.

The `run_simulations.sh` script combines the phase-defining and the
stage-defining mdp files to create the run file used by `gmx grompp`.

#### GROMACS topology files
The topology files for the different systems mentioned are stored in
folders named by their PDB ID in `systems/`. The base topologies
(named `system.top`) do not contain the intramolecular interaction
which restrains the water to the protein, because these restraints
vary based on the stage in the thermodynamic cycle and the simulation
phase. These interactions are stored in separate files in the same
directories named `waterRestraintLambdaIndependent.top` and
`waterRestraintLambdaDependent.top`, respectively.

#### Configuration files
The input configurations (named `system.gro`) used for the
minimization are stored in the same folders as their respective
topology files.

## Tests
### `tests/test_run_simulations.sh`
This script tests the `scripts/run_simulations.sh` and
`scripts/prepare_nes_structures.sh` files by making sure that the
typical usage (including calls to GROMACS) runs without errors. The
test is invoked by calling `bash tests/test_run_simulations.sh` from
the root of the repository. The test script expects `gmx` to be in the
path.

A GitHub action runs this test on every push and every PR to `main`.

### Python tests
The Python files located in the `water_nes/` folder are covered by a
combination of unit and regression tests. They can be invoked by
```shell
python -m pytest tests
```

The tests require `pytest` and `pytest-regressions` to be installed.

The regression tests are storing the output of previous runs of the
tests, and are comparing the output of the current run to this
reference output. As is the nature of regression tests, this means
that some changes will (and should!) make these tests fail.
The regression tests allow us to monitor exactly what changes
happened, and check whether they were intended. If such a
change is intended, the reference data can be regenerated using
```shell
python -m pytest tests --force-regen
```

## HPC3 files
The scripts located at `hpc3/` are example files on how to use the scripts
in this repository on UC Irvine's HPC3 cluster. These files are mostly
stored to improve reproducibility within the Mobley Lab, but might also be
useful as templates to others.
