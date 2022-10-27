import argparse

import MDAnalysis as mda
import MDAnalysis.transformations as mda_transformations
import numpy as np


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


def find_support(
    reference,
    pocket_selection,
    init_residue,
    final_residue,
    min_dist_pocket,
    min_dist_others,
):
    pocket_cog = reference.select_atoms(pocket_selection).center_of_geometry()
    backbone = reference.select_atoms(f"backbone")

    for r1 in range(init_residue, final_residue):
        if backbone.residues[r1] in reference.select_atoms(pocket_selection).residues:
            continue
        p1 = backbone.residues[r1].atoms.select_atoms("backbone").positions[1, :]
        # Is it far enough from pocket?
        if np.linalg.norm(p1 - pocket_cog) < min_dist_pocket:
            continue
        for r2 in range(r1 + 1, final_residue):
            if (
                backbone.residues[r2]
                in reference.select_atoms(pocket_selection).residues
            ):
                continue
            p2 = backbone.residues[r2].atoms.select_atoms("backbone").positions[1, :]
            # Is it far enough from pocket and p1?
            if np.linalg.norm(p2 - pocket_cog) < min_dist_pocket:
                continue
            if np.linalg.norm(p1 - p2) < min_dist_others:
                continue
            for r3 in range(r2 + 1, final_residue):
                if (
                    backbone.residues[r3]
                    in reference.select_atoms(pocket_selection).residues
                ):
                    continue
                p3 = (
                    backbone.residues[r3].atoms.select_atoms("backbone").positions[1, :]
                )
                # Is it far enough from pocket and p1 / p2?
                if np.linalg.norm(p3 - pocket_cog) < min_dist_pocket:
                    continue
                if np.linalg.norm(p1 - p3) < min_dist_others:
                    continue
                if np.linalg.norm(p2 - p3) < min_dist_others:
                    continue
                distances = np.linalg.norm(
                    [
                        p1 - pocket_cog,
                        p2 - pocket_cog,
                        p3 - pocket_cog,
                        p1 - p2,
                        p1 - p3,
                        p2 - p3,
                    ],
                    axis=1,
                )

                a1 = backbone.residues[r1].atoms.select_atoms("backbone").ids[1] + 1
                a2 = backbone.residues[r2].atoms.select_atoms("backbone").ids[1] + 1
                a3 = backbone.residues[r3].atoms.select_atoms("backbone").ids[1] + 1

                best_tuple = (a1, a2, a3)
                print(best_tuple)
                print(distances)
                return best_tuple


def print_posre(file_object, pocket_ids, support_atoms, file_type, force_constant):
    print("[ position_restraints ]", file=file_object)
    if file_type == "OFF_ON":
        print(
            f";{'idx':>4} {'fun':>3} {'fcx_A':>10} {'fcy_A':>10} {'fcz_A':>10}"
            f" {'fcx_B':>10} {'fcy_B':>10} {'fcz_B':>10}",
            file=file_object,
        )
    elif file_type == "ON" or file_type == "OFF":
        print(
            f";{'idx':>4} {'fun':>3} {'fcx':>10} {'fcy':>10} {'fcz':>10}",
            file=file_object,
        )
    else:
        raise NotImplemented(f"Unknown file type {file_type}")

    func = 1
    fcx = force_constant
    fcy = force_constant
    fcz = force_constant

    if file_type == "OFF_ON":
        for idx in pocket_ids:
            print(
                f"{idx:>5d} {func:>3d} {0:>10.2f} {0:>10.2f} {0:>10.2f} {fcx:>10.2f} {fcy:>10.2f} {fcz:>10.2f}",
                file=file_object,
            )
        for idx in support_atoms:
            print(
                f"{idx:>5d} {func:>3d} {fcx:>10.2f} {fcy:>10.2f} {fcz:>10.2f} {fcx:>10.2f} {fcy:>10.2f} {fcz:>10.2f}",
                file=file_object,
            )
    elif file_type == "ON":
        for idx in pocket_ids:
            print(
                f"{idx:>5d} {func:>3d} {fcx:>10.2f} {fcy:>10.2f} {fcz:>10.2f}",
                file=file_object,
            )
        for idx in support_atoms:
            print(
                f"{idx:>5d} {func:>3d} {fcx:>10.2f} {fcy:>10.2f} {fcz:>10.2f}",
                file=file_object,
            )
    elif file_type == "OFF":
        for idx in support_atoms:
            print(
                f"{idx:>5d} {func:>3d} {fcx:>10.2f} {fcy:>10.2f} {fcz:>10.2f}",
                file=file_object,
            )
    else:
        raise NotImplementedError(f"Unknown file type {file_type}")


def transform_restraint_units(value, in_units, out_units):
    transformation = {
        "kcal/mol/Å^2": {
            "kcal/mol/Å^2": 1.0,
            "kJ/mol/nm^2": 4.184 * 100,
        },
        "kJ/mol/nm^2": {
            "kcal/mol/Å^2": 0.239006 * 1e-2,
            "kJ/mol/nm^2": 1.0,
        },
    }

    return value * transformation[in_units][out_units]


def write_pocket_restraint_topology(
    tpr_file: str,
    gro_file: str,
    output_directory: str,
) -> None:
    reference = load_universe(tpr_file, gro_file, transfer_to_memory=True)
    pocket_selection = get_pocket_selection_string(reference)
    pocket_ids = reference.select_atoms(pocket_selection).ids + 1

    backbone = reference.select_atoms(f"backbone")

    min_dist_pocket = 15
    min_dist_others = 7.5

    # support_atoms = find_support(
    #     reference,
    #     pocket_selection,
    #     4,
    #     len(backbone.residues) - 4,
    #     min_dist_pocket,
    #     min_dist_others,
    # )

    support_atoms = reference.select_atoms("protein and name CA").ids + 1

    for file_type in ["ON", "OFF", "OFF_ON"]:
        with open(
            f"{output_directory}/posre_pocket_{file_type.lower()}.itp", "w"
        ) as out_file:
            print_posre(
                file_object=out_file,
                pocket_ids=pocket_ids,
                support_atoms=support_atoms,
                file_type=file_type,
                force_constant=transform_restraint_units(
                    0.5, "kcal/mol/Å^2", "kJ/mol/nm^2"
                ),
            )

    with open(f"{output_directory}/posre_protein.itp", "w") as out_file:
        print_posre(
            file_object=out_file,
            pocket_ids=reference.select_atoms("protein and (not name H*)").ids + 1,
            support_atoms=[],
            file_type="ON",
            force_constant=1000,
        )

    with open(f"{output_directory}/posre_ligand.itp", "w") as out_file:
        ligand_ids = reference.select_atoms(
            "resname LIG and (not name H*) and (not name DH*)"
        ).ids
        ligand_ids -= np.min(ligand_ids)
        ligand_ids += 1
        print_posre(
            file_object=out_file,
            pocket_ids=ligand_ids,
            support_atoms=[],
            file_type="ON",
            force_constant=1000,
        )


def command_line_entry_point():
    parser = argparse.ArgumentParser(description="Write pocket restraint files.")
    parser.add_argument("tpr_file", type=str, help="The minimization tpr file")
    parser.add_argument(
        "input_configuration", type=str, help="The minimized configuration"
    )
    parser.add_argument(
        "output_directory",
        type=str,
        help="A directory to write the restraint files to",
    )

    args = parser.parse_args()
    write_pocket_restraint_topology(
        tpr_file=args.tpr_file,
        gro_file=args.input_configuration,
        output_directory=args.output_directory,
    )


if __name__ == "__main__":
    command_line_entry_point()
