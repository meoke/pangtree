"""Parameters of the Affinity Tree generation algorithms."""

from io import StringIO
from pathlib import Path
from typing import Union


class Blosum(object):
    """BLOSUM matrix file content and its full path.

    This file is used by the poa software to determine consensuses paths
    in poagraph. The matrix in this file must contain symbol used by
    pangtreebuild for missing nucleotides/proteins. Lower-case is interpreted
    as nucleotides and upper-case as proteins.

    Args:
        file_content: Blosum file content required to validate its content.
        filepath: Path to the Blosum file, required by the poa software.

    Raises:
        ValueError: If file is empty or the BLOSUM matrix cannot be found.

    Attributes:
        filepath (Path): Full path to the blosum file.
            Required by poa software.
    """

    def __init__(self, file_content: StringIO, filepath: Path):
        self.filepath = filepath
        self._blosum_lines = file_content.readlines()
        self._blosum_symbols_line = None
        self._read_and_raise_exception_if_incorrect()

    def _read_and_raise_exception_if_incorrect(self) -> None:
        """Blosum file validation.

        Raises:
            ValueError: If file is empty or the BLOSUM matrix cannot be found.
        """

        if not self._blosum_lines:
            raise ValueError("""Empty blosum file. Provide a valid blosum file
                                or use te default one.""")

        blosum_symbols_line = None
        for blosum_line in self._blosum_lines:
            if len(blosum_line) > 0 and blosum_line[0] == " ":
                blosum_symbols_line = blosum_line
                break
        if not blosum_symbols_line:
            raise ValueError("""Cannot find the horizontal line of symbols
                                in blosum file. It should be the first line
                                beginning with a whitespace.""")
        self._blosum_symbols_line = blosum_symbols_line

    def check_if_symbol_is_present(self, missing_symbol: str) -> bool:
        """Checks if missing_base is present in the matrix.

        Args:
            missing_symbol: The symbol to be searched for in the BLOSUM matrix.

        Returns:
            bool: True if the symbol is included in the BLOSUM matrix.
                Otherwise an exception is raised.

        Raises:
            ValueError: If missing_base has length other than 1 or the symbol
                is not included in the BLOSUM matrix.
        """

        blosum_symbols = str.strip(self._blosum_symbols_line).split(" ")
        if len(missing_symbol) != 1:
            raise ValueError(f"""The missing symbol must be single character.
                                 Got{len(missing_symbol)} instead.""")
        if missing_symbol in blosum_symbols:
            return True
        else:
            raise ValueError("""Cannot find symbol used as missing base
                                in blosum file.""")


class Hbmin(object):
    """(Poa parameter) Minimum value of sequence compatibility to consensus.

    Args:
        value: Hbmin value. Must be in range [0,1].

    Raises:
        ValueError: If hbmin value is not in [0,1].
        TypeError: If value is not floatable.

    Attributes:
        value (float): Hbmin value.
    """

    def __init__(self, value: Union[str, float] = None):
        self.value: float = float(value) if value is not None else 0.9
        self._raise_exception_if_incorrect(self.value)

    def _raise_exception_if_incorrect(self, value: float) -> None:
        """Check if Hbmin value is correct, raise Exception if not.

        Args:
            value: Hbmin value to be checked.
        Raises:
            TypeError: If value is not float or not convertable to float.
            ValueError: If value is not in range [0,1].
        """

        try:
            v = float(value)
        except ValueError:
            raise TypeError(f"{value} was passed, a float excpected.")
        if v < 0 or v > 1:
            raise ValueError(f"Hbmin must be in range [0,1].")


class Stop(object):
    """Minimum value of AffinityTree leaves mincomp.

    Args:
        value: Stop value.

    Raises:
        ValueError: If Stop value is not in range [0,1].

    Attributes:
        value (float): The stop value.
    """

    def __init__(self, value: Union[str, float] = None):
        self.value: float = float(value) if value is not None else 0.99
        self._raise_exception_if_incorrect(self.value)

    def _raise_exception_if_incorrect(self, value: float) -> None:
        if value < 0:
            raise ValueError("STOP must be greater than 0.")
        if value > 1:
            raise ValueError("STOP must be smaller than 1.")


class P(object):
    """Parameter that changes compatiblity linear meaning to compatibility**P.

    Args:
        value: P value.

    Attributes:
        value (float): P value.
    """

    def __init__(self, value: Union[str, float] = None):
        self.value: float = float(value) if value is not None else 1
