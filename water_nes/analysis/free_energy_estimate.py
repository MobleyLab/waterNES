from dataclasses import asdict, dataclass


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
