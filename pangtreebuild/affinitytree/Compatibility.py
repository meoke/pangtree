"""Compatibility storage."""
from typing import Union, Any

from pangtreebuild.affinitytree import parameters


class Compatibility(object):
    """Asymetric similiarity measure of two poagraph paths.

    Args:
        compatibility: Raw compatibility value - count of common nodes devided by length of one of the paths.
        p: Parameter to control compatibility value interpretation. Compatibility is raised to the power of P.

    Attributes:
        value (float): Compatibility value raised to the power of p.
        p (float): P parameter value.

    """

    def __init__(self, compatibility: float, p: parameters.P = parameters.P(1)):
        self.value: float = compatibility**p.value
        self.p: float = p.value

    def _check_p_equality(self, other: Union["Compatibility", Any]) -> None:
        if isinstance(other, Compatibility):
            assert self.p == other.p, 'Cannot compare compatibilities with different p values.'
        else:
            return

    def __eq__(self, other: "Compatibility") -> bool:
        self._check_p_equality(other)
        return self.value == other.value

    def __lt__(self, other: "Compatibility") -> bool:
        self._check_p_equality(other)
        return self.value < other.value

    def __le__(self, other: "Compatibility") -> bool:
        self._check_p_equality(other)
        return self.value <= other.value

    def __gt__(self, other: "Compatibility") -> bool:
        self._check_p_equality(other)
        return self.value > other.value

    def __ge__(self, other: "Compatibility") -> bool:
        self._check_p_equality(other)
        return self.value >= other.value

    def __sub__(self, other: "Compatibility") -> "Compatibility":
        self._check_p_equality(other)
        return Compatibility(self.value - other.value, parameters.P(self.p))

    def __str__(self) -> str:
        return f"""{self.value}"""

    def __repr__(self) -> str:
        return f"""value: {self.value}, p: {self.p}"""

    def base_value(self) -> "Compatibility":
        """Get compatibility value without P transformation.

        Returns:
        Compatibility object with the original compatibility value and P=1."""

        return Compatibility(self.value ** (1 / self.p))
