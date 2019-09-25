from pangtreebuild.affinitytree import parameters


class Compatibility:
    """Asymetric similiarity measure of two poagraph paths.

    Attributes:
        value (float): compatibility value raised to the power of p
        p (float): P parameter value

    """

    def __init__(self, compatibility: float, p: parameters.P = parameters.P(1)):
        """
        Creates Compatibility object.

        Arguments:
        compatibility: Raw compatibility value - count of common nodes devided by length of one of the paths.
        p: Parameter to control compatibility value interpretation. Compatibility is raised to the power of P.
        """

        self.value: float = compatibility**p.value
        self.p: float = p.value

    def _check_p_equality(self, other):
        assert self.p == other.p, 'Cannot compare compatibilities with different p values.'

    def __eq__(self, other):
        self._check_p_equality(other)
        return self.value == other.value

    def __lt__(self, other):
        self._check_p_equality(other)
        return self.value < other.value

    def __le__(self, other):
        self._check_p_equality(other)
        return self.value <= other.value

    def __gt__(self, other):
        self._check_p_equality(other)
        return self.value > other.value

    def __ge__(self, other):
        self._check_p_equality(other)
        return self.value >= other.value

    def __sub__(self, other):
        self._check_p_equality(other)
        return Compatibility(self.value - other.value, parameters.P(self.p))

    def __str__(self):
        return f"""{self.value}"""

    def __repr__(self):
        return f"""value: {self.value}, p: {self.p}"""

    def base_value(self):
        """Get compatibility value without P transformation.

        Returns:
        Compatibility object with the original compatibility value."""

        return Compatibility(self.value ** (1 / self.p))