from io import StringIO
from typing import Optional, Tuple, List
from pathlib import Path


class ConsensusInputError(Exception):
    pass


class Blosum:
    """File with BLOSUM matrix.
    This file is used by poa software to determine consensuses paths in poagraph.
    The matrix in this file must contain symbol for missing nucleotides/proteins.
    Lower-case is interpreted as nucleteotides and upper-case as proteins."""

    def __init__(self, file_content: StringIO, filepath: Optional[Path], missing_base_symbol: str):
        self.filepath = filepath
        self._raise_exception_if_incorrect(file_content, missing_base_symbol)
        self.filecontent = file_content

    @staticmethod
    def _raise_exception_if_incorrect(filecontent: StringIO, missing_base_symbol: str):
        blosum_lines = filecontent.readlines()
        if not blosum_lines:
            raise ConsensusInputError("Empty blosum file. Provide a valid blosum file or use te default one.")

        blosum_symbols_line = None
        for blosum_line in blosum_lines:
            if len(blosum_line) > 0 and blosum_line[0] == " ":
                blosum_symbols_line = blosum_line
                break
        if not blosum_symbols_line:
            raise ConsensusInputError("Cannot find the horizontal line of symbols in blosum file. "
                                      "It should begin with a whitespace.")
        blosum_symbols = str.strip(blosum_symbols_line).split(" ")

        if missing_base_symbol in blosum_symbols:
            return True
        else:
            raise ConsensusInputError("Cannot find symbol used as missing base in blosum file.")


class Hbmin:
    """The minimum value of sequence compatibility to generated consensus."""

    def __init__(self, value: str):
        self._raise_exception_if_incorrect(value)
        self.value: float = float(value)

    def _raise_exception_if_incorrect(self, value: str) -> None:
        try:
            v = float(value)
        except ValueError:
            raise ConsensusInputError(f"{value} was passed, a float excpected.")
        if v < 0 or v > 1:
            raise ConsensusInputError(f"Hbmin must be in range [0,1].")


class Range:
    """Specify what part of sorted capabilities should be searched for node cutoff. E.g. [0.2,0.8]"""

    def __init__(self, value: List[str]):
        self.value: Tuple[float, float] = (float(value[0]), float(value[1]))
        self._raise_exception_if_incorrect(self.value)

    def _raise_exception_if_incorrect(self, value: Tuple[float, float]) -> None:
        if len(value) != 2:
            raise Exception("CUTOFF SEARCH RANGE must have length 2.")
        if value[1] < value[0]:
            raise Exception("CUTOFF SEARCH RANGE first value must be smaller or equal to second value.")
        if value[0] < 0 \
                or value[0] > 1 \
                or value[1] < 0 \
                or value[1] > 1:
            raise Exception("CUTOFF SEARCH RANGE values must be in the range of [0,1].")

