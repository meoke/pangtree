import abc
from enum import Enum
from io import StringIO

from pangtreebuild.datamodel.Node import Base
from pangtreebuild.datamodel.Sequence import SequenceID


class FastaProviderException(Exception):
    pass


class FastaProviderOption(Enum):
    """todo"""

    NO = 0
    NCBI = 1
    FILE = 2


class UseCache:
    """Used if Fasta Provider is \"NCBI\"
       Sequences downloaded from NCBI are stored and
       reused between uses of this program."""
    pass


class FastaProvider(abc.ABC):
    """Poagraph building from maf file parameter. Ignored when -not_dag parameter is set.
       Maf file usually contains not full sequences but only parts of them, aligned to each other.
       To build an exact poagraph the full sequences must be retrieved from: ncbi or local file system."""

    @abc.abstractmethod
    def get_base(self, sequence_id: SequenceID, i: int) -> Base:
        pass

    @staticmethod
    def get_raw_sequence_from_fasta(fasta_handle: StringIO) -> str:
        _ = fasta_handle.readline()
        return fasta_handle.read().replace('\n', '')
