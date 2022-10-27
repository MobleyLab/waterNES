import argparse

import MDAnalysis as mda
import MDAnalysis.transformations as mda_transformations
import numpy as np
import scipy.optimize


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
    distance = 5
    pocket_resids = universe.select_atoms(
        f"backbone and around {distance} group trapped_water_o",
        trapped_water_o=trapped_water_o,
    ).residues.resids

    # Increase distance if we don't have at least 3 residues
    while len(pocket_resids) < 3:
        distance += 1
        pocket_resids = universe.select_atoms(
            f"backbone and around {distance} group trapped_water_o",
            trapped_water_o=trapped_water_o,
        ).residues.resids

    # Use only non-H atoms for selection
    resids_selection = " or ".join([f"resid {resid}" for resid in pocket_resids])
    pocket_selection = f"({resids_selection}) and (not name H*)"
    return pocket_selection


def distance_to_virtual_site(positions1, positions2, target_positions, parameter):
    virtual_site = (1 - parameter) * positions1 + parameter * positions2
    return np.linalg.norm(target_positions - virtual_site, axis=1)


def difference_with_virtual_site(positions1, positions2, target_positions, parameter):
    virtual_site = (1 - parameter) * positions1 + parameter * positions2
    return target_positions - virtual_site


def find_virtual_site(
    pocket_positions,
    target_positions,
):
    min_res = None
    min_tuple = (None, None)
    min_parameter = None

    num_pocket_atoms = pocket_positions.shape[1]

    for atom_idx_1 in range(num_pocket_atoms):
        for atom_idx_2 in range(atom_idx_1 + 1, num_pocket_atoms):

            def target_function(parameter_a):
                return distance_to_virtual_site(
                    pocket_positions[:, atom_idx_1, :],
                    pocket_positions[:, atom_idx_2, :],
                    target_positions,
                    parameter_a,
                ).sum()

            res = scipy.optimize.minimize_scalar(
                target_function, bounds=(0, 1), method="Bounded"
            )
            if min_res is None or (res.success and res.fun < min_res):
                min_res = res.fun
                min_tuple = (atom_idx_1, atom_idx_2)
                min_parameter = res.x

    return min_tuple[0], min_tuple[1], min_parameter, min_res


def find_virtual_site_and_update_configuration(
    tpr_file: str,
    input_configuration: str,
    output_configuration: str,
    input_topology: str,
    output_topology: str,
) -> None:
    reference = load_universe(tpr_file, input_configuration, transfer_to_memory=True)

    pocket_selection = get_pocket_selection_string(reference)
    vs_atom_1, vs_atom_2, vs_parameter, _, = find_virtual_site(
        np.array([reference.select_atoms(pocket_selection).positions]),
        reference.select_atoms("resname MOL and name O").positions,
    )

    gromacs_atom1 = reference.select_atoms(pocket_selection)[vs_atom_1].id + 1
    gromacs_atom2 = reference.select_atoms(pocket_selection)[vs_atom_2].id + 1
    gromacs_virtual = reference.select_atoms("protein").ids[-1] + 2
    print("VIRTUAL SITE MIN:")

    print("[ virtual_sites2 ]")
    print(f";{'atom_nr':>8} {'atom_i':>8} {'atom_j':>8} {'fct':>8} {'a':>8}")
    print(
        f" {gromacs_virtual:>8} {gromacs_atom1:>8}"
        f" {gromacs_atom2:>8} {1:>8} {vs_parameter:>8.4f}"
    )

    virtual_site_position = (1 - vs_parameter) * reference.select_atoms(
        pocket_selection
    )[vs_atom_1].position.flatten() + vs_parameter * reference.select_atoms(
        pocket_selection
    )[
        vs_atom_2
    ].position.flatten()

    with open(input_configuration) as in_conf:
        with open(output_configuration, "w") as out_conf:
            virtual_site_written = False
            for num_line, line in enumerate(in_conf):
                if num_line == 1:
                    line = f"{int(line.strip()) + 1}\n"
                if not virtual_site_written and "LIG" in line:
                    out_conf.write(
                        f"{0:5}{'ATT':<3}{'AT1':>7}{gromacs_virtual:5} "
                        f"{virtual_site_position[0]:7.3f} "
                        f"{virtual_site_position[1]:7.3f} "
                        f"{virtual_site_position[2]:7.3f}\n"
                    )
                    virtual_site_written = True
                out_conf.write(line)

    # Note: This assumes that the protein is the first `moleculetype` block,
    #       and MOL the second
    with open(input_topology) as in_top:
        topology = []
        molecule_block_counter = 0
        for line in in_top:
            if "[ moleculetype ]" == line.strip():
                molecule_block_counter += 1
            if molecule_block_counter == 1 and "[ bonds ]" == line.strip():
                while topology[-1].strip() == "":
                    topology.pop(-1)
                last_residue = int(topology[-1].split()[2])
                topology.extend(
                    [
                        f"#ifdef SOL_RESTRAINT_ON\n",
                        f"{gromacs_virtual:5}  attachX {last_residue+1}   "
                        f"ATT   AT1  {gromacs_virtual:5}   "
                        f"0.00000   0.00000  ; dummy atom - no mass\n",
                        f"#endif\n",
                        f"#ifdef SOL_RESTRAINT_OFF\n",
                        f"{gromacs_virtual:5}  attach  {last_residue+1}   "
                        f"ATT   AT1  {gromacs_virtual:5}   "
                        f"0.00000   0.00000  ; dummy atom - no mass\n",
                        f"#endif\n",
                        f"#ifdef SOL_RESTRAINT_OFF_ON\n",
                        f"{gromacs_virtual:5}  attach  {last_residue+1}   "
                        f"ATT   AT1  {gromacs_virtual:5}   "
                        f"0.00000   0.00000  attachX  0.00000   0.00000\n",
                        f"#endif\n",
                        f"\n",
                    ]
                )
            topology.append(line)
        # TODO: Put virtual site at end of block

    topology.extend(
        [
            f"\n",
            f"#ifdef MOL_RESTRAINT_ON\n",
            f"[ intermolecular_interactions ]\n",
            f"[ bonds ]\n",
            f";      ai       aj     type       bA       kA       bB       kB\n",
            f"    {gromacs_virtual:5}    {gromacs_virtual+1:5}        "
            f"6      0.0    209.2      0.0    209.2\n",
            f"#endif\n",
            f"#ifdef MOL_RESTRAINT_OFF_ON\n",
            f"[ intermolecular_interactions ]\n",
            f"[ bonds ]\n",
            f";      ai       aj     type       bA       kA       bB       kB\n",
            f"    {gromacs_virtual:5}    {gromacs_virtual+1:5}        "
            f"6      0.0      0.0      0.0    209.2\n",
            f"#endif\n",
        ]
    )

    with open(output_topology, "w") as out_top:
        for line in topology:
            out_top.write(line)


def command_line_entry_point():
    parser = argparse.ArgumentParser(
        description="Find definition of virtual site based on minimized configuration."
    )
    parser.add_argument("tpr_file", type=str, help="The minimization tpr file")
    parser.add_argument(
        "input_configuration", type=str, help="The minimized configuration"
    )
    parser.add_argument(
        "output_configuration",
        type=str,
        help="A file name to write the minimized configuration "
        "including the virtual site to",
    )
    parser.add_argument(
        "input_topology",
        type=str,
        help="The input topology",
    )
    parser.add_argument(
        "output_topology",
        type=str,
        help="A file name to write the altered topology to",
    )

    args = parser.parse_args()
    find_virtual_site_and_update_configuration(
        tpr_file=args.tpr_file,
        input_configuration=args.input_configuration,
        output_configuration=args.output_configuration,
        input_topology=args.input_topology,
        output_topology=args.output_topology,
    )


if __name__ == "__main__":
    command_line_entry_point()
