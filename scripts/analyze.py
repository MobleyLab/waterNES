import argparse
import os
import pickle
import shutil
import sys
import warnings
from typing import Dict, Iterable, List, Optional, Tuple

import alchemlyb
import MDAnalysis as mda
import numpy as np
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.postprocessors.units import get_unit_converter
from alchemlyb.preprocessing import slicing, statistical_inefficiency
from MDAnalysis import transformations as mda_transformations
from MDAnalysis.analysis import distances as mda_distances
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader
from pymbar.timeseries import ParameterError

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


def calculate_mbar_inner(
    stages: Iterable,
    cycle_directory: str,
    drop_first: bool,
    file_modifier: Optional[str],
    subsample: bool,
    conservative: bool = True,
) -> MBAR:
    u_nk_all_stages = []
    file_name = f"dhdl_{file_modifier}.xvg" if file_modifier is not None else "dhdl.xvg"
    for counter, stage in enumerate(stages):
        with warnings.catch_warnings():
            # Tons of pandas warnings here
            warnings.simplefilter(action="ignore", category=FutureWarning)
            full_u_nk = extract_u_nk(
                f"{cycle_directory}/stage{stage}/prod/{file_name}", T=298.15
            )

            if drop_first:
                full_u_nk = full_u_nk.drop((0.0, 0.0, 0.0), axis=1)

            u_nk = slicing(full_u_nk, lower=2000)
            if subsample:
                u_nk_all_stages.append(
                    statistical_inefficiency(
                        u_nk,
                        series=u_nk[u_nk.columns[counter]],
                        conservative=conservative,
                    )
                )
            else:
                u_nk_all_stages.append(u_nk)

    return MBAR().fit(alchemlyb.concat(u_nk_all_stages))


def calculate_mbar(
    stages: Iterable,
    cycle_directory: str,
    drop_first: bool,
    file_modifier: Optional[str],
) -> MBAR:
    """
    TODO: Need to better understand this. MBAR sometimes fails when subsampling,
          and sometimes when _not_ subsampling. The results look reasonable either
          way, so for now I work around this by rerunning the analysis without
          subsampling if needed. This shouldn't affect the current analysis, but it
          would be good to understand the reason for it and possibly fix it!
    """
    try:
        return calculate_mbar_inner(
            stages,
            cycle_directory,
            drop_first,
            file_modifier,
            subsample=True,
            conservative=False,
        )
    except ParameterError:
        print("Non-conservative subsampling failed, trying conservative subsampling")
    try:
        return calculate_mbar_inner(
            stages,
            cycle_directory,
            drop_first,
            file_modifier,
            subsample=True,
            conservative=True,
        )
    except ParameterError:
        print("Conservative subsampling failed, trying without subsampling")

    return calculate_mbar_inner(
        stages, cycle_directory, drop_first, file_modifier, subsample=False
    )


def get_conversion_factor(input_units: str, output_units: str) -> float:
    if input_units != "kJ/mol/nm^2":
        raise NotImplementedError(f"Unknown input units {input_units}")
    if output_units == "kJ/mol/nm^2":
        return 1
    elif output_units == "kcal/mol/nm^2":
        return 0.239006
    elif output_units == "kJ/mol/Å^2":
        return 0.01
    elif output_units == "kcal/mol/Å^2":
        return 0.239006 * 0.01

    raise NotImplementedError(f"Unknown output units {output_units}")


def get_pocket_restraint_force_constant(cycle_directory: str, units: str) -> float:
    restraint_force_constant = None
    with open(f"{cycle_directory}/../posre_pocket_on.itp") as itp_file:
        for line in itp_file:
            if "position_restraints" in line or line.startswith(";"):
                continue
            restraint_force_constant = float(line.split()[-1])  # kJ/mol/nm^2
            break
    assert restraint_force_constant is not None
    return restraint_force_constant * get_conversion_factor(
        input_units="kJ/mol/nm^2", output_units=units
    )


def get_water_restraint_force_constant(cycle_directory: str, units: str) -> float:
    restraint_force_constant = None
    with open(f"{cycle_directory}/../system.top") as top_file:
        block_found = False
        for line in top_file:
            if "intermolecular_interactions" in line:
                block_found = True
                continue
            if block_found and not line.startswith(";") and len(line.split()) == 7:
                restraint_force_constant = float(line.split()[-1])  # kJ/mol/nm^2
                break

    assert restraint_force_constant is not None
    return restraint_force_constant * get_conversion_factor(
        input_units="kJ/mol/nm^2", output_units=units
    )


def calculate_lower_edge_free_energy(
    cycle_directory: str, file_modifier: Optional[str]
) -> Dict[str, FreeEnergyEstimate]:
    output_units = "kcal/mol"
    if file_modifier is None:
        stages = ["7", "6", "6.1", "6.2", "6.3", "6.4", "6.5", "6.6", "6.7", "5"]
    else:
        stages = ["7", "6", "5"]

    mbar = calculate_mbar(
        stages, cycle_directory, drop_first=True, file_modifier=file_modifier
    )
    edge_f = sum(
        [
            get_unit_converter(output_units)(mbar.delta_f_).iloc[idx][idx - 1]
            for idx in range(2, mbar.delta_f_.shape[0])
        ]
    )
    edge_f_error = sum(
        [
            get_unit_converter(output_units)(mbar.d_delta_f_).iloc[idx][idx - 1]
            for idx in range(2, mbar.d_delta_f_.shape[0])
        ]
    )

    restraint_force_constant = get_water_restraint_force_constant(
        cycle_directory=cycle_directory, units="kcal/mol/nm^2"
    )

    rt = 1.9872159e-3 * 298.15  # kcal/mol/K * K
    edge_h = rt * np.log(
        55.5 * 0.6022 * (2 * np.pi * rt / restraint_force_constant) ** (3 / 2)
    )
    # Unit analysis:
    # (2 * np.pi * rt / restraint_force_constant) ** (3 / 2) --> nm^3
    # Molar concentration of water:
    #   55.5 mol / L == 55.5 mol / (10^24 nm^3) == 55.5 * 6.022 * 10^23 / (10^24 nm^3)
    #                == 55.5 * 0.6022 nm^-3

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


def calculate_upper_edge_free_energy(
    cycle_directory: str, file_modifier: Optional[str]
) -> Dict[str, FreeEnergyEstimate]:
    output_units = "kcal/mol"
    if file_modifier is None:
        stages = ["1", "2", "3", "3.1", "3.2", "3.3", "3.4", "3.5", "3.6", "3.7", "4"]
    else:
        stages = ["1", "2", "3", "4"]

    mbar = calculate_mbar(
        stages, cycle_directory, drop_first=False, file_modifier=file_modifier
    )
    edge_d = sum(
        [
            get_unit_converter(output_units)(mbar.delta_f_).iloc[idx - 1][idx]
            for idx in range(3, mbar.delta_f_.shape[0])
        ]
    )
    edge_d_error = sum(
        [
            get_unit_converter(output_units)(mbar.d_delta_f_).iloc[idx - 1][idx]
            for idx in range(3, mbar.d_delta_f_.shape[0])
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


def calculate_nes_edges(cycle_directory: str) -> Dict[str, FreeEnergyEstimate]:
    edge_and_stages = [
        ("E", 4, 5),
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

        try:
            (
                free_energies[f"Edge {edge}"],
                free_energies[f"Work edge {edge}"],
            ) = calculate_nes_free_energy(
                xvg_files_forward_transition=xvg_files_forward,
                xvg_files_backward_transition=xvg_files_backward,
                temperature=298.15,
                output_units="kcal/mol",
                bootstrapping_repeats=0,
            )
        except ValueError:
            free_energies[f"Edge {edge}"] = FreeEnergyEstimate(
                value=np.nan, error=np.nan, units="kcal/mol", bootstrap_error=0
            )
            free_energies[f"Work edge {edge}"] = None

    return free_energies


def simple_restraint(force_constant: float, distance: np.ndarray) -> np.array:
    distance = np.array(distance)
    return 0.5 * force_constant * (distance ** 2)


def read_xvg(xvg_file_name: str) -> Tuple[Dict[str, np.ndarray], List[str]]:
    contents = {}
    header = []
    with open(xvg_file_name, "r") as xvg_file:
        for line in xvg_file:
            if line.startswith("#") or line.startswith("@"):
                header.append(line)
                continue
            entries = line.split()
            time = entries[0]
            assert time not in contents
            contents[time] = np.array([float(entry) for entry in entries[1:]])
    num_frames = len(contents)
    return contents, header


def write_xvg(
    input_xvg_file: str,
    output_xvg_file: str,
    water_restraint_energy: Optional[np.ndarray],
    position_restraint_energy: Optional[np.ndarray],
    current_lambda_state: int,
    alternative_xvg_file: Optional[str] = None,
) -> None:
    # Read existing xvg file
    input_xvg_contents, header = read_xvg(xvg_file_name=input_xvg_file)
    # This is a bit brittle (i.e., depending on our exact simulation setup),
    # but for now, we'll just hardcode the xvg file format
    entry = {
        "potential": 0,
        "vdw-lambda": 1,
        "bonded-lambda": 2,
        "restraint-lambda": 3,
        "foreign lambda 0": 4,
        "foreign lambda 1": 5,
        "foreign lambda 2": 6,
        "foreign lambda 3": 7,
        "foreign lambda 4": 8,
        "foreign lambda 5": 9,
        "foreign lambda 6": 10,
        "foreign lambda 7": 11,
        "foreign lambda 8": 12,
        "foreign lambda 9": 13,
        "foreign lambda 10": 14,
        "pV": 15,
    }
    num_entries = len(input_xvg_contents)
    if water_restraint_energy is not None:
        assert len(water_restraint_energy) == num_entries
        assert current_lambda_state == 0
    if position_restraint_energy is not None:
        assert len(position_restraint_energy) == num_entries

    if position_restraint_energy is not None:
        header.remove(
            '@ s5 legend "\\xD\\f{}H \\xl\\f{} to (0.0000, 1.0000, 0.0000)"\n'
        )
        header.remove(
            '@ s6 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.0000)"\n'
        )
        header.remove(
            '@ s7 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.0100)"\n'
        )
        header.remove(
            '@ s8 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.0200)"\n'
        )
        header.remove(
            '@ s9 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.0400)"\n'
        )
        header.remove(
            '@ s10 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.0800)"\n'
        )
        header.remove(
            '@ s11 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.1600)"\n'
        )
        header.remove(
            '@ s12 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.3200)"\n'
        )
        header.remove(
            '@ s13 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 0.6400)"\n'
        )
        header.remove(
            '@ s14 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 1.0000)"\n'
        )
        header.remove('@ s15 legend "pV (kJ/mol)"\n')
        header.append(
            '@ s5 legend "\\xD\\f{}H \\xl\\f{} to (1.0000, 1.0000, 1.0000)"\n'
        )
        header.append('@ s6 legend "pV (kJ/mol)"\n')

    if (alternative_xvg_file is not None) and (current_lambda_state == 10):
        # The alternative file leaves stages 4 and 5 unchanged, removes only the
        # intermediate states.
        with open(alternative_xvg_file, "w") as out_xvg:
            out_xvg.write("".join(line for line in header))
            for idx, time in enumerate(input_xvg_contents):
                xvg_line = input_xvg_contents[time]
                for lambda_idx in range(1, 10):
                    xvg_line[entry[f"foreign lambda {lambda_idx}"]] = np.nan
                out_xvg.write(
                    f"{time} "
                    + " ".join(
                        f"{entry:#.8g}" for entry in xvg_line if not np.isnan(entry)
                    )
                    + "\n"
                )

    with open(output_xvg_file, "w") as out_xvg:
        out_xvg.write("".join(line for line in header))

        for idx, time in enumerate(input_xvg_contents):
            xvg_line = input_xvg_contents[time]
            for lambda_idx in range(1, 10):
                xvg_line[entry[f"foreign lambda {lambda_idx}"]] = np.nan

            if water_restraint_energy is not None:
                old_value = xvg_line[entry["bonded-lambda"]]
                new_value = water_restraint_energy[idx]
                xvg_line[entry["bonded-lambda"]] = new_value
                for lambda_idx in [10]:
                    xvg_line[entry[f"foreign lambda {lambda_idx}"]] += (
                        new_value - old_value
                    )

            if position_restraint_energy is not None:
                for lambda_idx in range(3, 10):
                    xvg_line[entry[f"foreign lambda {lambda_idx}"]] = np.nan
                old_value = xvg_line[entry["restraint-lambda"]]
                new_value = position_restraint_energy[idx]
                xvg_line[entry["restraint-lambda"]] = new_value
                if current_lambda_state == 10:
                    xvg_line[entry["potential"]] += new_value - old_value
                    for lambda_idx in [0]:
                        xvg_line[entry[f"foreign lambda {lambda_idx}"]] += (
                            old_value - new_value
                        )
                else:
                    xvg_line[entry[f"foreign lambda 10"]] += new_value - old_value

            out_xvg.write(
                f"{time} "
                + " ".join(f"{entry:#.8g}" for entry in xvg_line if not np.isnan(entry))
                + "\n"
            )

    if (alternative_xvg_file is not None) and (current_lambda_state != 10):
        # The alternative file is identical to the new file for other lambda states
        shutil.copy2(output_xvg_file, alternative_xvg_file)


def do_analysis(
    system: str,
    cycle_directory: str,
    read_file: bool,
    do_distances: bool,
    postprocess_trapped: bool,
    postprocess_pocket: bool,
    do_free_energy: bool,
):
    reference = load_universe(
        topology=f"{cycle_directory}/stage1/min/topol.tpr",
        trajectory=f"{cycle_directory}/../minimized.gro",
    )
    pocket_selection_string = get_pocket_selection_string(reference)

    analysis = {}
    if read_file:
        with open(f"analysis{system}.pickle", "rb") as in_file:
            analysis = pickle.load(in_file)

    if postprocess_trapped and "distances" not in analysis:
        do_distances = True

    if do_distances:
        distances = {}
        all_stages = [1, 4, 5, 6, 7]
        for stage in all_stages:
            trajectory = load_universe(
                topology=f"{cycle_directory}/stage{stage}/prod/topol.tpr",
                trajectory=f"{cycle_directory}/stage{stage}/prod/traj_comp.xtc",
            )
            distances[f"stage{stage}"] = calculate_distances(
                universe=trajectory, pocket_selection_string=pocket_selection_string
            )
        analysis["distances"] = distances

    trapped_restraint_energy = None
    if postprocess_trapped:
        restraint_force_constant = get_water_restraint_force_constant(
            cycle_directory=cycle_directory, units="kJ/mol/Å^2"
        )
        trapped_restraint_energy = simple_restraint(
            force_constant=restraint_force_constant,
            distance=analysis["distances"]["stage1"]["closest water"],
        )

    position_restraint_energy = None
    if postprocess_pocket:
        all_stages = [1, 4, 5, 6, 7]
        position_restraint_energy = {}
        for stage in all_stages:
            # Load trajectory
            trajectory = load_universe(
                topology=f"{cycle_directory}/stage{stage}/prod/topol.tpr",
                trajectory=f"{cycle_directory}/stage{stage}/prod/traj_comp.xtc",
                transfer_to_memory=False,
            )
            # Align frames, and calculate positional restraint
            position_restraint_energy[f"stage{stage}"] = np.zeros(
                len(trajectory.trajectory)
            )
            force_constant = get_pocket_restraint_force_constant(
                cycle_directory=cycle_directory, units="kJ/mol/Å^2"
            )
            for (
                frame_idx,
                _,
            ) in enumerate(trajectory.trajectory):
                mda.analysis.align.alignto(
                    trajectory.select_atoms(pocket_selection_string),
                    reference.select_atoms(pocket_selection_string),
                )
                distances = np.linalg.norm(
                    trajectory.select_atoms(pocket_selection_string).positions
                    - reference.select_atoms(pocket_selection_string).positions,
                    axis=1,
                )
                position_restraint_energy[f"stage{stage}"][
                    frame_idx
                ] = simple_restraint(
                    force_constant=force_constant,
                    distance=distances,
                ).sum()

    if postprocess_trapped or postprocess_pocket:
        if postprocess_pocket:
            all_stages = [1, 4, 5, 6, 7]
        else:
            all_stages = [1]
        for stage in all_stages:
            # Write modified xvg file
            position_restraint_energy_for_stage = None
            if postprocess_pocket:
                position_restraint_energy_for_stage = position_restraint_energy[
                    f"stage{stage}"
                ]
            stage_to_lambda_state = {
                "stage1": 0,
                "stage4": 10,
                "stage5": 10,
                "stage6": 2,
                "stage7": 1,
            }
            write_xvg(
                input_xvg_file=f"{cycle_directory}/stage{stage}/prod/dhdl.xvg",
                output_xvg_file=f"{cycle_directory}/stage{stage}/prod/dhdl_mod.xvg",
                water_restraint_energy=trapped_restraint_energy if stage == 1 else None,
                position_restraint_energy=position_restraint_energy_for_stage,
                current_lambda_state=stage_to_lambda_state[f"stage{stage}"],
                alternative_xvg_file=f"{cycle_directory}/stage{stage}/prod/dhdl_alt.xvg"
                if stage_to_lambda_state[f"stage{stage}"] == 10
                else None,
            )

    if do_free_energy:
        free_energies = {
            "Edge I": FreeEnergyEstimate(value=0, error=0, units="kcal/mol"),
            "Edge J": FreeEnergyEstimate(value=-6.13, error=0.01, units="kcal/mol"),
        }
        free_energies_modified = {
            "Edge I": FreeEnergyEstimate(value=0, error=0, units="kcal/mol"),
            "Edge J": FreeEnergyEstimate(value=-6.13, error=0.01, units="kcal/mol"),
        }
        free_energies_alternative = {
            "Edge I": FreeEnergyEstimate(value=0, error=0, units="kcal/mol"),
            "Edge J": FreeEnergyEstimate(value=-6.13, error=0.01, units="kcal/mol"),
        }
        free_energies.update(
            calculate_upper_edge_free_energy(cycle_directory, file_modifier=None)
        )
        free_energies.update(
            calculate_lower_edge_free_energy(cycle_directory, file_modifier=None)
        )
        free_energies_modified.update(
            calculate_upper_edge_free_energy(cycle_directory, file_modifier="mod")
        )
        free_energies_modified.update(
            calculate_lower_edge_free_energy(cycle_directory, file_modifier="mod")
        )
        free_energies_alternative.update(
            calculate_upper_edge_free_energy(cycle_directory, file_modifier="alt")
        )
        free_energies_alternative.update(
            calculate_lower_edge_free_energy(cycle_directory, file_modifier="alt")
        )
        nes_edges = calculate_nes_edges(cycle_directory)
        free_energies.update(nes_edges)
        free_energies_modified.update(nes_edges)
        free_energies_alternative.update(nes_edges)
        analysis["free energy"] = free_energies
        analysis["free energy modified"] = free_energies_modified
        analysis["free energy alternative"] = free_energies_alternative

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
    parser.add_argument(
        "--postProcessTrapped",
        action="store_true",
        help="Postprocess the trapped water restraint in stage 1",
    )
    parser.add_argument(
        "--postProcessPocket",
        action="store_true",
        help="Postprocess the binding pocket restraints",
    )
    parser.add_argument("--fe", action="store_true", help="Do free energy calculations")
    args = parser.parse_args()
    do_analysis(
        system=args.system,
        cycle_directory=args.directory,
        read_file=args.read,
        do_distances=args.dist,
        postprocess_trapped=args.postProcessTrapped,
        postprocess_pocket=args.postProcessPocket,
        do_free_energy=args.fe,
    )


if __name__ == "__main__":
    command_line_entry_point()
