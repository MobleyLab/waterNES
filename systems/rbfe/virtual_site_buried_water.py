#!/usr/bin/env python
# coding: utf-8

# This is just to make sure we start in the root of the repo - for this script we only need to make sure that we can access the `systems/` folder with the different ligands.

# In[5]:


import os
import sys

os.chdir("..")
sys.path.append(os.getcwd())
print(os.getcwd())  # Should be root of repo!


# This is a bunch of code to 
# * read in MDAnalyis universes
# * define a selection which returns backbone atoms that represent the binding pocket of the water molecule, defined by proximity to the trapped water
# * find a virtual site by looping over a set of atoms (the binding pocket defined above), and finding the pair of atoms which allows to place a virtual site closest to the target (the trapped water)
# 
# Note that the exact results of this will depend on how the pocket is defined. By default, the pocket is defined as any backbone atom which is part of a residue which has a at least one backbone atom witin 5A of the target. If that selection does not yield at least 3 residues, increase the radius until at least 3 residues are part of the selection. Choosing these parameters differently would yield different virtual sites. But the definition of the virtual site has proven to be somewhat flexible, as long as a _light_ restraining force constant is chosen.

# In[7]:


import numpy as np
import scipy.optimize
import MDAnalysis as mda
import MDAnalysis.transformations as mda_transformations


def add_workflow(universe):
    protein_and_ligand = universe.select_atoms("protein or resname LIG")
    water_and_virtual_site = universe.select_atoms("resname MOL or resname ATT")

    center_selection = universe.select_atoms("resname LIG")
    if len(center_selection) == 0:
        center_selection = universe.select_atoms("resname MOL")

    assert len(center_selection) > 0

    workflow = [
        mda_transformations.unwrap(universe.atoms),
        mda_transformations.center_in_box(center_selection, center="mass"),
        mda_transformations.wrap(protein_and_ligand, compound="fragments"),
        mda_transformations.wrap(water_and_virtual_site, compound="atoms"),
    ]
    universe.trajectory.add_transformations(*workflow)

    return universe


def load_universe(topology, trajectory, transfer_to_memory=True, step=1, workflow=True):
    universe = mda.Universe(topology, trajectory)
    if transfer_to_memory:
        universe.transfer_to_memory(step=step)
    if workflow:
        return add_workflow(universe)
    return universe


def get_pocket_selection_string(universe):
    # Find binding pocket selection
    trapped_water_o = universe.select_atoms("resname MOL and name O")
    #print(trapped_water_o)
#    trapped_water_o = universe.select_atoms("resname MOL")
    distance = 8
    pocket_resids = universe.select_atoms(
        f"backbone and around {distance} group trapped_water_o",
        trapped_water_o=trapped_water_o,
    ).residues.resids
    #print(pocket_resids)
    #breakpoint()
    # Increase distance if we don't have at least 3 residues
    while len(pocket_resids) < 3:
        distance += 1
        pocket_resids = universe.select_atoms(
            f"backbone and around {distance} group trapped_water_o",
            trapped_water_o=trapped_water_o,
        ).residues.resids

    # Use only non-H atoms for selection
    resids_selection = " or ".join([f"resid {resid}" for resid in pocket_resids])
    pocket_selection = f"({resids_selection}) and backbone"
    #breakpoint()
    return pocket_selection


def distance_to_virtual_site(positions1, positions2, target_positions, parameter):
    virtual_site = (1 - parameter) * positions1 + parameter * positions2
    return np.linalg.norm(target_positions - virtual_site, axis=1)


def find_virtual_site(
    pocket_positions,
    target_positions,
):
    min_res = None
    min_tuple = (None, None)
    min_parameter = None

    num_pocket_atoms = pocket_positions.shape[0]

    for atom_idx_1 in range(num_pocket_atoms):
        for atom_idx_2 in range(atom_idx_1 + 1, num_pocket_atoms):

            def target_function(parameter_a):
                return distance_to_virtual_site(
                    pocket_positions[atom_idx_1, :],
                    pocket_positions[atom_idx_2, :],
                    target_positions,
                    parameter_a,
                ).sum()

            res = scipy.optimize.minimize_scalar(
                target_function, bounds=(0, 1), method="Bounded"
            )
            if res.success and (min_res is None or res.fun < min_res):
                min_res = res.fun
                min_tuple = (atom_idx_1, atom_idx_2)
                min_parameter = res.x

    return min_tuple[0], min_tuple[1], min_parameter, min_res


# To test whether the above definition is indeed the definition used in Pascal's simulations, compare the calculated virtual sites to the ones found in the topology.

# In[3]:


# Check all ABFE systems:
proteins = {
    'S8-S11' : ['S8-S11']
}
abfe_ligands = [
    ligand
    for protein in proteins
    for ligand in proteins[protein]
]
abfe_ligands


# In[4]:


for ligand in abfe_ligands:
    tpr_file = sys.argv[1].strip()
    input_configuration = sys.argv[2].strip()
    topology = sys.argv[3].strip()
    reference = load_universe(tpr_file, input_configuration, transfer_to_memory=True)
    pocket = reference.select_atoms(get_pocket_selection_string(reference))
    trapped_water_o = reference.select_atoms("resname MOL and name O")

    # grep -A2 "virtual_sites2" topology
    print("Virtual site in system.top:")
    with open(topology) as read_topology:
        for line in read_topology:
            if "virtual_sites2" in line:
                print(f"  {line}", end="")
                print(f"  {next(read_topology)}", end="")
                topology_vsite = next(read_topology)
                print(f"  {topology_vsite}", end="")
                break
    ai = int(topology_vsite.split()[1])
    aj = int(topology_vsite.split()[2])
    param = float(topology_vsite.split()[4])
    distance = distance_to_virtual_site(
        reference.select_atoms('all').positions[ai-1],
        reference.select_atoms('all').positions[aj-1],
        trapped_water_o.positions,
        param
    )[0]
    print(f"  Distance: {distance:.2f}nm\n")
    
    result = find_virtual_site(
        pocket.positions,
        trapped_water_o.positions
    )
    # Note: MDAnalysis atoms are 0-based, GROMACS atoms are 1-based
    print("Calculated virtual site")
    print(f"  atom_i = {pocket.ids[result[0]] + 1} | atom_j = {pocket.ids[result[1]] + 1} | "
          f"a = {result[2]:.4f} | Distance = {result[3]:.2f}nm\n")


# We can see that the calculated virtual sites perfectly agree with the ones in the topology with the exception of 4STD and 7STD. It's likely that I used a slightly bigger cutoff to define the pocket for these cases - maybe 8A instead of 5A. Otherwise, the above code is indeed defining the virtual sites used in my simulations.
