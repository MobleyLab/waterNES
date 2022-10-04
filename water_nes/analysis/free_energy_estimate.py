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

    def __str__(self):
        if self.bootstrap_error != 0.0:
            error = self.bootstrap_error
        else:
            error = self.error

        return f"{self.value:7.3f} +- {error:5.3f} {self.units}"

    def __round__(self, num_digits: int = 0) -> "FreeEnergyEstimate":
        return FreeEnergyEstimate(
            value=round(self.value, num_digits),
            error=round(self.error, num_digits),
            units=self.units,
            bootstrap_error=round(self.bootstrap_error, num_digits),
        )

    def as_dict(self) -> dict:
        return asdict(self)

    def __neg__(self) -> "FreeEnergyEstimate":
        return FreeEnergyEstimate(
            value=-self.value,
            error=self.error,
            units=self.units,
            bootstrap_error=self.bootstrap_error,
        )

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
