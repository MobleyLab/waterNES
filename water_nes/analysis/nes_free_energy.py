from typing import Any, Dict, List, Tuple, Union

import numpy as np
from pmx.analysis import read_dgdl_files
from pmx.estimators import BAR

from .free_energy_estimate import FreeEnergyEstimate


def make_list_of_lists(
    list_or_list_of_lists_1: Union[List[List[Any]], List[Any]],
    list_or_list_of_lists_2: Union[List[List[Any]], List[Any]],
) -> Tuple[List[List[Any]], List[List[Any]]]:
    r"""Unify input to single type

    :param list_or_list_of_lists_1: A list of items, or a list of lists of items
    :param list_or_list_of_lists_2: A list of items, or a list of lists of items
    :return: Both input parameters as list of lists of items
    """
    if not isinstance(list_or_list_of_lists_1, List) or not isinstance(
        list_or_list_of_lists_2, List
    ):
        raise TypeError("Unsupported type")

    # Put input and results in lists for easier traversal
    xvg_files_types = [None, None]
    xvg_files = [list_or_list_of_lists_1, list_or_list_of_lists_2]

    # Test items of both lists
    for idx in range(2):
        if len(xvg_files[idx]) > 0 and all(
            isinstance(item, List) for item in xvg_files[idx]
        ):
            # All items of the list are lists
            xvg_files_types[idx] = "List of lists"
        elif not any(isinstance(item, List) for item in xvg_files[idx]):
            # None of the items of the list are a list
            xvg_files_types[idx] = "Simple list"

    # Return if both lists have the same, compatible type
    if xvg_files_types[0] == xvg_files_types[1]:
        if xvg_files_types[0] == "Simple list":
            return [list_or_list_of_lists_1], [list_or_list_of_lists_2]
        if xvg_files_types[0] == "List of lists":
            return list_or_list_of_lists_1, list_or_list_of_lists_2

    raise TypeError("Unsupported type")


def find_unique_list_size(dictionary: Dict[str, List[List[Any]]]) -> int:
    r"""Return the size of the innermost list

    :param dictionary: A dictionary containing list of lists
    :return: The unique size of the inner lists, if it exists
    """
    size = None
    for key in dictionary:
        for inner_list in dictionary[key]:
            if size is None:
                size = len(inner_list)
            elif size != len(inner_list):
                raise RuntimeError("No common inner size list")
    return size


def convert_energy(energy_in_kj_mol: float, units: str) -> float:
    r"""Convert energy value from kJ/mol to different units

    :param energy_in_kj_mol: Energy value in kJ/mol
    :param units: Unit to convert to, currently supports 'kcal/mol' and 'kJ/mol'
    :return: Energy value in new units
    """
    if units == "kJ/mol":
        return energy_in_kj_mol
    if units == "kcal/mol":
        return energy_in_kj_mol * 0.2390057361
    raise NotImplementedError(f"Unit {units} not implemented.")


def calculate_nes_free_energy(
    xvg_files_forward_transition: Union[List[List[Any]], List[Any]],
    xvg_files_backward_transition: Union[List[List[Any]], List[Any]],
    temperature: float,
    output_units: str,
    bootstrapping_repeats: int = 0,
) -> FreeEnergyEstimate:
    r"""Calculate free energy estimate from swarm of non-equilibrium switching simulations

    :param xvg_files_forward_transition:
        A list of xvg files (single transition), or a list of lists of xvg files
        (multiple transitions) for stage A
    :param xvg_files_backward_transition:
        A list of xvg files (single transition), or a list of lists of xvg files
        (multiple transitions) for stage B
    :param temperature:
        The temperature at which the simulations were performed, in Kelvin
    :param output_units:
        The units to return the resulting free energy, either 'kcal/mol' or 'kJ/mol'
    :param bootstrapping_repeats:
        The number of bootstrapping repeats to perform for error estimating
        (default: 0, no error estimate)
    :return: The calculated free energy
    """

    # Unify input to be of type List[List[Any]] to simplify further use
    xvg_files = {}
    try:
        (xvg_files["forward"], xvg_files["backward"]) = make_list_of_lists(
            xvg_files_forward_transition, xvg_files_backward_transition
        )
    except TypeError:
        raise TypeError(
            "At least one of xvg_files_forward_transition or "
            "xvg_files_backward_transition was neither a list nor a list of lists, "
            "or the two arguments did not have the same type"
        )

    # Find size of simulation swarm
    try:
        swarm_size = find_unique_list_size(xvg_files)
    except RuntimeError:
        raise TypeError(
            "All lists in xvg_files_forward_transition and "
            "xvg_files_backward_transition need to have the same length."
        )

    # This assumes that the forward transition starts at lambda=0,
    # and the backward transition at lambda=1
    # Could be made an input parameter if the need arises
    start_lambda = {"forward": 0, "backward": 1}

    # For each simulation in the swarm, sum up work of the different sub transitions
    # The result is an array of size `swarm_size` for each direction with the work
    # of the full transition
    work = {"forward": np.zeros(swarm_size), "backward": np.zeros(swarm_size)}
    for direction in work:
        for files_per_transition in xvg_files[direction]:
            # read_dgdl_files return type is wrongly detected, so supress the warning
            # noinspection PyTypeChecker
            work[direction] += np.array(
                read_dgdl_files(
                    files_per_transition,
                    lambda0=start_lambda[direction],
                    invert_values=False,
                )
            )

    # Use all the work estimates to make a BAR Estimate
    estimate = BAR(
        work["forward"], work["backward"], T=temperature, nboots=bootstrapping_repeats
    )

    return FreeEnergyEstimate(
        value=convert_energy(estimate.dg, output_units),
        error=convert_energy(estimate.err, output_units),
        bootstrap_error=convert_energy(estimate.err_boot, output_units)
        if bootstrapping_repeats > 0
        else 0,
        units=output_units,
    )
