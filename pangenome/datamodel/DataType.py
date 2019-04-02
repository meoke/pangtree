from enum import Enum
from typing import NewType, Dict

from .Sequence import SequenceID, Sequence


class DataType(Enum):
    """todo"""

    Nucleotides = 0
    Proteins = 1


Sequences = NewType("Sequences", Dict[SequenceID, Sequence])
