#!/bin/bash
#SBATCH --account=DMOBLEY_LAB
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --constraint=avx512

# The source command can be replaced by loading a module (if the GROMACS version you
# want to use is available) or by a pointer to any other GROAMCS installation.
# The next line is telling shellcheck not to worry about the sourced file.
# shellcheck source=/dev/null
source ~/bin/gmx2022.2+/bin/GMXRC

gmx mdrun -ntomp 1 -cpi state.cpt
