from io import StringIO
from typing import Optional, Tuple, List, Union
from pathlib import Path

from pangtreebuild.datamodel.input_types import MissingSymbol


class ConsensusInputError(Exception):
    pass


class Blosum:
    """File with BLOSUM matrix.
    This file is used by poa software to determine consensuses paths in poagraph.
    The matrix in this file must contain symbol for missing nucleotides/proteins.
    Lower-case is interpreted as nucleteotides and upper-case as proteins."""

    def __init__(self, file_content: StringIO, filepath: Path):
        self.filepath = filepath
        self.blosum_lines = file_content.readlines()
        self.blosum_symbols_line = None
        self._raise_exception_if_incorrect()

    def _raise_exception_if_incorrect(self):
        if not self.blosum_lines:
            raise ConsensusInputError("Empty blosum file. Provide a valid blosum file or use te default one.")

        blosum_symbols_line = None
        for blosum_line in self.blosum_lines:
            if len(blosum_line) > 0 and blosum_line[0] == " ":
                blosum_symbols_line = blosum_line
                break
        if not blosum_symbols_line:
            raise ConsensusInputError("Cannot find the horizontal line of symbols in blosum file. "
                                      "It should be the first line beginning with a whitespace.")
        self.blosum_symbols_line = blosum_symbols_line

    def check_if_symbol_is_present(self, missing_symbol: str):
        """Checks if missng__symbol is present in the matrix."""

        blosum_symbols = str.strip(self.blosum_symbols_line).split(" ")
        if missing_symbol in blosum_symbols:
            return True
        else:
            raise ConsensusInputError("Cannot find symbol used as missing base in blosum file.")


class Hbmin:
    """The minimum value of sequence compatibility to generated consensus."""

    def __init__(self, value: Union[str, float] = None):
        self.value: float = float(value) if value is not None else 0.9
        self._raise_exception_if_incorrect(self.value)

    def _raise_exception_if_incorrect(self, value: float) -> None:
        try:
            v = float(value)
        except ValueError:
            raise ConsensusInputError(f"{value} was passed, a float excpected.")
        if v < 0 or v > 1:
            raise ConsensusInputError(f"Hbmin must be in range [0,1].")


class Stop:
    """Value of node compatibility above which the node is no more split."""
    
    def __init__(self, value: Union[str, float] = None):
        self.value: float = float(value) if value is not None else 0.99
        self._raise_exception_if_incorrect(self.value)

    def _raise_exception_if_incorrect(self, value: float) -> None:
        if value < 0:
            raise ConsensusInputError("STOP must be greater than 0.")
        if value > 1:
            raise ConsensusInputError("STOP must be smaller than 1.")


class P:
    """Value that changes compatiblity linear meaning to compatibility**P."""
    """Any floatable value is allowed."""

    def __init__(self, value: Union[str, float] = None):
        self.value: float = float(value) if value is not None else 1



