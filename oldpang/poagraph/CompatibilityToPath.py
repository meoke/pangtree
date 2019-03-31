class CompatibilityToPath:
    def __init__(self, compatibility: float, p: float):
        self.value: float = compatibility**p

    def __eq__(self, other):
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __le__(self, other):
        return self.value <= other.value

    def __gt__(self, other):
        return self.value > other.value

    def __ge__(self, other):
        return self.value >= other.value

    def __sub__(self, other):
        return CompatibilityToPath(self.value - other.value, 1)

    def __str__(self):
        return f"""{self.value}"""

    def __repr__(self):
        return f"""{self.value}"""
