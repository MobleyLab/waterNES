import argparse
import numpy as np
from MDAnalysis.analysis import distances as mda_distances

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


def vdw_repulsive_energy(c12: float, distance: np.ndarray) -> np.ndarray:
    distance = np.array(distance)
    return c12 * (distance ** (-12))


def get_solvent_restraint_c12(cycle_directory: str, units: str) -> float:
    restraint_c12 = None
    with open(f"{cycle_directory}/../system.top") as top_file:
        block_found = False
        for line in top_file:
            if "nonbond_params" in line:
                block_found = True
                continue
            if block_found and not line.startswith(";") and len(line.split()) == 5:
                fields = line.split()
                assert fields[0] == "attachX"
                sigma = float(fields[-2])
                epsilon = float(fields[-1])
                assert sigma < 0  # This signifies a purely repulsive potential in gmx
                restraint_c12 = 4 * epsilon * (sigma ** 12)
                break

    assert restraint_c12 is not None
    return restraint_c12 * get_conversion_factor(
        input_units="kJ/mol*nm^12", output_units=units
    )


def simple_restraint(force_constant: float, distance: np.ndarray) -> np.ndarray:
    distance = np.array(distance)
    return 0.5 * force_constant * (distance ** 2)

def get_conversion_factor(input_units: str, output_units: str) -> float:
    conversions = {
        "kJ/mol/nm^2": {
            "kJ/mol/nm^2": 1,
            "kcal/mol/nm^2": 0.239006,
            "kJ/mol/Å^2": 0.01,
            "kcal/mol/Å^2": 0.239006 * 0.01,
        },
        "kJ/mol*nm^12": {
            "kJ/mol*nm^12": 1,
            "kcal/mol*nm^12": 0.239006,
            "kJ/mol*Å^12": 1e12,
            "kcal/mol*Å^12": 0.239006 * 1e12,
        },
    }

    if input_units not in conversions:
        raise NotImplementedError(f"Unknown input units {input_units}")
    if output_units not in conversions[input_units]:
        raise NotImplementedError(
            f"Unknown output units {output_units} " f"for input unit {input_units}"
        )

    return conversions[input_units][output_units]


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
    attachment = universe.select_atoms("resname ATT")
    trapped_water = universe.select_atoms("resname MOL and name O")
    solvent = universe.select_atoms("resname HOH  and name O")
#    all_waters = universe.select_atoms("(resname HOH or resname MOL) and name O")
    pocket = universe.select_atoms(pocket_selection_string)
    distance_dict = {
        "trapped water": [],
        "solvent water 1st": [],
        "solvent water 2nd": [],
        "closest water": [],
        "closest solvent": [],
        "second closest solvent": [],
        "third closest solvent": [],
        "fourth closest solvent": [],
        "fifth closest solvent": [],
        "pocket": [],
    }
    for _ in universe.trajectory:
        distance_trapped_water = mda_distances.distance_array(
            attachment.positions,
            trapped_water.positions,
            box=universe.dimensions,
            backend='serial'
        ).flatten()
        distance_solvent = np.sort(
            mda_distances.distance_array(
                attachment.positions,
                solvent.positions,
                box=universe.dimensions,
                backend='serial'
            ).flatten()
        )
#        distance_all_waters= np.sort(
#            mda_distances.distance_array(
#                attachment.positions,
#                all_waters.positions,
#                box=universe.dimensions,
#            ).flatten()
#        )

        distance_pocket = mda_distances.distance_array(
            attachment.positions,
            pocket.positions,
            box=universe.dimensions,
        )

        distance_dict["trapped water"].append(distance_trapped_water)
        distance_dict["solvent water 1st"].append(distance_solvent[0])
        distance_dict["solvent water 2nd"].append(distance_solvent[1])

        distance_solvent = distance_solvent[:5]
        distance_solvent = np.sort(np.concatenate((distance_solvent, distance_trapped_water), axis=0))


#        distance_dict["closest water"].append(
#            min(distance_trapped_water, distance_solvent[0])
#        )
        distance_dict["closest water"].append(distance_solvent[0])
        distance_dict["closest solvent"].append(distance_solvent[1])
        distance_dict["second closest solvent"].append(distance_solvent[2])
        distance_dict["third closest solvent"].append(distance_solvent[3])
        distance_dict["fourth closest solvent"].append(distance_solvent[4])
        distance_dict["fifth closest solvent"].append(distance_solvent[5])
#        distance_dict["second closest solvent"].append(np.array2string(distance_solvent[2:], precision=2, separator=','))
#        print(type(np.array2string(distance_solvent[2:], precision=2, separator=',')), distance_dict["second closest solvent"])
        distance_dict["pocket"].append(distance_pocket)

    for key in distance_dict:
        distance_dict[key] = np.array(distance_dict[key])

    return distance_dict


def command_line_entry_point():
    parser = argparse.ArgumentParser(description=
            "This script analyzes the thermodynamic cycle for the relative "
            "binding free energy calculation of buried water systems."
            "The script was originally written by Pascal Merz (DE Shaw research),"
            "and was further edited by Swapnil Wagle (University of California Irvine, swapnilw@uci.edu)"
            "To run the script, either type: 'python analysis.py -h', or type: ")

    parser.add_argument("system", type=str, help="The name of the system. Output files will be generated using the system name")
    parser.add_argument("directory", type=str, help=(
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
        help="Postprocess the trapped water and solvent restraint in stage 1",
    )
    parser.add_argument(
        "--postProcessPocket",
        action="store_true",
        help="Postprocess the binding pocket restraints",
    )
    parser.add_argument(
        "--fe_option",
        type=str,
        help="Do free energy calculations for the stages. There are the following four options to do that "
        "1: Stages: 1, 4, 5, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6, 7"
        "2: Stages: 1, 2, 3, 4, 5, 6, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6, 7" 
        "3: Stages: 1, 2, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 4, 5, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6, 7"
        "4: Custom option: Stages should be given in the format"
    )
    return parser.parse_args()

