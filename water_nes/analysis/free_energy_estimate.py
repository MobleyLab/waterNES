from typing import NamedTuple


class FreeEnergyEstimate(NamedTuple):
    r"""
    Represents a free energy estimate, including information
    on error estimates and units
    """
    value: float
    error: float
    bootstrap_error: float
    units: str

    def __round__(self, num_digits: int = 0) -> "FreeEnergyEstimate":
        return FreeEnergyEstimate(
            value=round(self.value, num_digits),
            error=round(self.error, num_digits),
            bootstrap_error=round(self.bootstrap_error, num_digits),
            units=self.units,
        )
