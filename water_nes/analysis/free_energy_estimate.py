from dataclasses import asdict, dataclass
from typing import List

import numpy as np


@dataclass
class FreeEnergyEstimate:
    r"""
    Represents a free energy estimate, including information
    on error estimates and units
    """
    value: float
    error: float
    units: str
    bootstrap_error: float = 0.0

    def __round__(self, num_digits: int = 0) -> "FreeEnergyEstimate":
        return FreeEnergyEstimate(
            value=round(self.value, num_digits),
            error=round(self.error, num_digits),
            units=self.units,
            bootstrap_error=round(self.bootstrap_error, num_digits),
        )

    def as_dict(self) -> dict:
        return asdict(self)

    def __add__(self, other: "FreeEnergyEstimate") -> "FreeEnergyEstimate":
        if self.units != other.units:
            raise NotImplementedError(
                "FreeEnergyEstimate objects with different units " "cannot be added"
            )
        return FreeEnergyEstimate(
            value=self.value + other.value,
            error=self.error + other.error,
            units=self.units,
            bootstrap_error=self.bootstrap_error + other.bootstrap_error,
        )

    def __sub__(self, other: "FreeEnergyEstimate") -> "FreeEnergyEstimate":
        if self.units != other.units:
            raise NotImplementedError(
                "FreeEnergyEstimate objects with different units "
                "cannot be subtracted"
            )
        return FreeEnergyEstimate(
            value=self.value - other.value,
            error=self.error + other.error,
            units=self.units,
            bootstrap_error=self.bootstrap_error + other.bootstrap_error,
        )

    def __mul__(self, other: float) -> "FreeEnergyEstimate":
        return FreeEnergyEstimate(
            value=self.value * other,
            error=self.error * other,
            units=self.units,
            bootstrap_error=self.bootstrap_error * other,
        )

    def __truediv__(self, other: float) -> "FreeEnergyEstimate":
        return FreeEnergyEstimate(
            value=self.value / other,
            error=self.error / other,
            units=self.units,
            bootstrap_error=self.bootstrap_error / other,
        )


def mean_free_energy(
    free_energy_estimates: List[FreeEnergyEstimate],
) -> FreeEnergyEstimate:
    r"""
    Returns a free energy estimate with value equal to the mean of a list of
    free energy estimates, and error equal to the standard deviation of the
    list of values.

    Parameters
    ----------
    free_energy_estimates : List[FreeEnergyEstimate]
        List of independent free energy estimates

    Returns
    -------
    FreeEnergyEstimate
        Mean of the list of free energy estimates, with its
        standard deviation as error

    """
    if not free_energy_estimates:
        raise TypeError("Cannot calculate the average of an empty list.")
    if any(
        [
            estimate.units != free_energy_estimates[0].units
            for estimate in free_energy_estimates
        ]
    ):
        raise NotImplementedError(
            "Averaging free energy estimates with different "
            "units is not implemented."
        )

    values = np.array([estimate.value for estimate in free_energy_estimates])
    return FreeEnergyEstimate(
        value=values.mean(),
        error=values.std(),
        units=free_energy_estimates[0].units,
    )
