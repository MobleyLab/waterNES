import argparse
import os
import pickle
import sys
import warnings

import alchemlyb
import MDAnalysis as mda
import numpy as np
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.postprocessors.units import get_unit_converter
from alchemlyb.preprocessing import slicing, statistical_inefficiency
from MDAnalysis import transformations as mda_transformations
from MDAnalysis.analysis import distances as mda_distances

sys.path.append(os.getcwd())
from water_nes.analysis.free_energy_estimate import FreeEnergyEstimate
from water_nes.analysis.nes_free_energy import calculate_nes_free_energy


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


def calculate_distances(universe, pocket_selection_string):
    reference = universe.select_atoms("resname ATT")
    trapped_water = universe.select_atoms("resname MOL and name O")
    solvent = universe.select_atoms("resname HOH and name O")
    pocket = universe.select_atoms(pocket_selection_string)
    distance_dict = {
        "trapped water": [],
        "solvent water 1st": [],
        "solvent water 2nd": [],
        "closest water": [],
        "pocket": [],
    }
    for _ in universe.trajectory:
        distance_trapped_water = mda_distances.distance_array(
            reference.positions,
            trapped_water.positions,
            box=universe.dimensions,
        ).flatten()
        distance_solvent = np.sort(
            mda_distances.distance_array(
                reference.positions,
                solvent.positions,
                box=universe.dimensions,
            ).flatten()
        )
        distance_pocket = mda_distances.distance_array(
            reference.positions,
            pocket.positions,
            box=universe.dimensions,
        )

        distance_dict["trapped water"].append(distance_trapped_water[0])
        distance_dict["solvent water 1st"].append(distance_solvent[0])
        distance_dict["solvent water 2nd"].append(distance_solvent[1])
        distance_dict["closest water"].append(
            min(distance_trapped_water[0], distance_solvent[0])
        )
        distance_dict["pocket"].append(distance_pocket)

    for key in distance_dict:
        distance_dict[key] = np.array(distance_dict[key])

    return distance_dict


def calculate_lower_edge_free_energy(cycle_directory):
    output_units = "kcal/mol"
    stages = ["7", "6", "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7", "5"]

    u_nk_all_stages = []
    for counter, stage in enumerate(stages):
        with warnings.catch_warnings():
            # Tons of pandas warnings here
            warnings.simplefilter(action='ignore', category=FutureWarning)
            u_nk = slicing(
                extract_u_nk(
                    f"{cycle_directory}/stage{stage}/prod/dhdl.xvg", T=298.15
                ).drop((0.0, 0.0, 0.0), axis=1),
                lower=2000,
            )
            u_nk_all_stages.append(
                statistical_inefficiency(u_nk, series=u_nk[u_nk.columns[counter]])
            )

    mbar = MBAR().fit(alchemlyb.concat(u_nk_all_stages))
    edge_f = sum(
        [
            get_unit_converter(output_units)(mbar.delta_f_).iloc[idx][idx - 1]
            for idx in range(2, len(u_nk_all_stages))
        ]
    )
    edge_f_error = sum(
        [
            get_unit_converter(output_units)(mbar.d_delta_f_).iloc[idx][idx - 1]
            for idx in range(3, len(u_nk_all_stages))
        ]
    )

    restraint_force_constant = None
    with open(f"{cycle_directory}/../system.top") as top_file:
        block_found = False
        for line in top_file:
            if "intermolecular_interactions" in line:
                block_found = True
                continue
            if block_found and not line.startswith(";") and len(line.split()) == 7:
                restraint_force_constant = (
                    float(line.split()[-1]) * 0.239006
                )  # kcal/mol/nm^2
    assert restraint_force_constant is not None

    rt = 1.9872159e-3 * 298.15  # kcal/mol/K * K
    edge_h = rt * np.log(
        55 * 0.6022 * (2 * np.pi * rt / restraint_force_constant) ** (3 / 2)
    )

    return {
        "Edge F": FreeEnergyEstimate(
            value=edge_f,
            error=edge_f_error,
            units=output_units,
        ),
        "Edge G": FreeEnergyEstimate(
            value=get_unit_converter(output_units)(mbar.delta_f_).iloc[1][0],
            error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[1][0],
            units=output_units,
        ),
        "Edge H": FreeEnergyEstimate(
            value=edge_h,
            error=0.0,
            units=output_units,
        ),
        "lower edge overlap": mbar.overlap_matrix,
    }


def calculate_upper_edge_free_energy(cycle_directory):
    output_units = "kcal/mol"
    stages = ["1", "2", "3", "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "4"]

    u_nk_all_stages = []
    for counter, stage in enumerate(stages):
        with warnings.catch_warnings():
            # Tons of pandas warnings here
            warnings.simplefilter(action='ignore', category=FutureWarning)
            u_nk = slicing(
                extract_u_nk(f"{cycle_directory}/stage{stage}/prod/dhdl.xvg", T=298.15),
                lower=2000,
            )
            u_nk_all_stages.append(
                statistical_inefficiency(u_nk, series=u_nk[u_nk.columns[counter]])
            )

    mbar = MBAR().fit(alchemlyb.concat(u_nk_all_stages))
    edge_d = sum(
        [
            get_unit_converter(output_units)(mbar.delta_f_).iloc[idx - 1][idx]
            for idx in range(3, len(u_nk_all_stages))
        ]
    )
    edge_d_error = sum(
        [
            get_unit_converter(output_units)(mbar.d_delta_f_).iloc[idx - 1][idx]
            for idx in range(3, len(u_nk_all_stages))
        ]
    )

    return {
        "Edge B": FreeEnergyEstimate(
            value=get_unit_converter(output_units)(mbar.delta_f_).iloc[0][1],
            error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[0][1],
            units=output_units,
        ),
        "Edge C": FreeEnergyEstimate(
            value=get_unit_converter(output_units)(mbar.delta_f_).iloc[1][2],
            error=get_unit_converter(output_units)(mbar.d_delta_f_).iloc[1][2],
            units=output_units,
        ),
        "Edge D": FreeEnergyEstimate(
            value=edge_d,
            error=edge_d_error,
            units=output_units,
        ),
        "upper edge overlap": mbar.overlap_matrix,
    }


def calculate_nes_edges(cycle_directory):
    edge_and_stages = [
        ("E", 4, 5),
        ("L", 3, 6),
        ("K", 2, 7),
    ]
    num_nes_repeats = 100

    free_energies = {}
    for (
        edge,
        stageA,
        stageB,
    ) in edge_and_stages:
        xvg_files_forward = [
            f"{cycle_directory}/stage{stageA}/NES/run{run}/dhdl.xvg"
            for run in range(1, num_nes_repeats + 1)
        ]
        xvg_files_backward = [
            f"{cycle_directory}/stage{stageB}/NES/run{run}/dhdl.xvg"
            for run in range(1, num_nes_repeats + 1)
        ]

        free_energies[f"Edge {edge}"] = calculate_nes_free_energy(
            xvg_files_forward_transition=xvg_files_forward,
            xvg_files_backward_transition=xvg_files_backward,
            temperature=298.15,
            output_units="kcal/mol",
            bootstrapping_repeats=0,
        )

    return free_energies


def do_analysis(system, cycle_directory, read_file, do_distances, do_free_energy):

    reference = load_universe(
        topology=f"{cycle_directory}/stage1/min/topol.tpr",
        trajectory=f"{cycle_directory}/../minimized.gro",
    )
    pocket_selection_string = get_pocket_selection_string(reference)

    analysis = {}
    if read_file:
        with open(f"analysis{system}.pickle", "rb") as in_file:
            analysis = pickle.load(in_file)

    if do_distances:
        distances = {}
        all_stages = [1, 2, 3, 4, 5, 6, 7]
        for stage in all_stages:
            trajectory = load_universe(
                topology=f"{cycle_directory}/stage{stage}/prod/topol.tpr",
                trajectory=f"{cycle_directory}/stage{stage}/prod/traj_comp.xtc",
            )
            distances[f"stage{stage}"] = calculate_distances(
                universe=trajectory, pocket_selection_string=pocket_selection_string
            )
        analysis["distances"] = distances

    if do_free_energy:
        free_energies = {
            "Edge I": FreeEnergyEstimate(value=0, error=0, units="kcal/mol"),
            "Edge J": FreeEnergyEstimate(value=-6.13, error=0.01, units="kcal/mol"),
        }
        free_energies.update(calculate_upper_edge_free_energy(cycle_directory))
        free_energies.update(calculate_lower_edge_free_energy(cycle_directory))
        free_energies.update(calculate_nes_edges(cycle_directory))
        analysis["free energy"] = free_energies

    with open(f"analysis{system}.pickle", "wb") as out_file:
        pickle.dump(analysis, out_file)


def command_line_entry_point():
    parser = argparse.ArgumentParser(description="Analyze ABFE for system")
    parser.add_argument("system", type=str, help="The name of the system")
    parser.add_argument(
        "directory",
        type=str,
        help=(
            "The top directory of the cycle to analyze, should contain folders "
            "named stageX where X is one of the valid stage numbers"
        ),
    )
    parser.add_argument(
        "--read",
        action="store_true",
        help="Read analysis{system}.pickle file (in current directory)",
    )
    parser.add_argument("--dist", action="store_true", help="Do distance calculations")
    parser.add_argument("--fe", action="store_true", help="Do free energy calculations")
    args = parser.parse_args()
    do_analysis(
        system=args.system,
        cycle_directory=args.directory,
        read_file=args.read,
        do_distances=args.dist,
        do_free_energy=args.fe,
    )


if __name__ == "__main__":
    command_line_entry_point()
