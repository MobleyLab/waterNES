from typing import Any

import alchemlyb
from alchemlyb.estimators import MBAR
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.postprocessors.units import get_unit_converter
from alchemlyb.preprocessing import slicing, statistical_inefficiency

from .free_energy_estimate import FreeEnergyEstimate


def calculate_endpoint_free_energy(
    file_lambda_0: Any,
    file_lambda_1: Any,
    start_time: int,
    end_time: int,
    output_units: str,
) -> FreeEnergyEstimate:
    r"""
    Calculate the free energy difference between two end points using (M)BAR

    Parameters
    ----------
    file_lambda_0 : str
        GROMACS xvg file at lambda=0
    file_lambda_1 : str
        GROMACS xvg file at lambda=1
    start_time : int
        Time at which to start considering values, in ps.
        Allows to discard equilibration time or to calculate dG for windows.
    end_time : int
        Time at which to stop considering values, in ps.
        Allows to calculate dG for windows.
    output_units : str
        Units to print and return the free energy estimate.
        One of 'kcal/mol', 'kJ/mol', 'kT'
    Returns
    -------
    FreeEnergyEstimate
        The estimate of the free energy difference between the two states.

    """
    # Read reduced potentials, and subsample them to reduce correlation
    # Note: Slicing before calling statistical_inefficiency is currently
    #       necessary: https://github.com/alchemistry/alchemlyb/issues/198
    u_nk_lambda_0 = slicing(
        extract_u_nk(file_lambda_0, T=298.15), lower=start_time, upper=end_time
    )
    u_nk_lambda_0_sub = statistical_inefficiency(
        u_nk_lambda_0, series=u_nk_lambda_0[u_nk_lambda_0.columns[0]]
    )

    u_nk_lambda_1 = slicing(
        extract_u_nk(file_lambda_1, T=298.15), lower=start_time, upper=end_time
    )
    u_nk_lambda_1_sub = statistical_inefficiency(
        u_nk_lambda_1, series=u_nk_lambda_1[u_nk_lambda_1.columns[-1]]
    )

    # Calculate MBAR
    mbar = MBAR().fit(alchemlyb.concat([u_nk_lambda_0_sub, u_nk_lambda_1_sub]))

    return FreeEnergyEstimate(
        value=get_unit_converter(output_units)(mbar.delta_f_).loc[0.0, 1.0].item(),
        error=get_unit_converter(output_units)(mbar.d_delta_f_).loc[0.0, 1.0].item(),
        bootstrap_error=0,
        units=output_units,
    )
