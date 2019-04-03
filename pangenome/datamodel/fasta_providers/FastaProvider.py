import abc
import re
from enum import Enum

from datamodel.Sequence import SequenceID


class FastaProviderException(Exception):
    pass


class FastaProviderOption(Enum):
    """todo"""

    NO = 0
    NCBI = 1
    FILE = 2


class EmailAddress:
    """E-mail address requiered when Fasta Provider Option is \"NCBI\"
       as Entrez API obligates the user to pass e-mail address."""

    def __init__(self, email_address: str):
        match = re.match(r'^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,4})$', email_address)
        if match is not None:
            self.value = email_address
        else:
            raise FastaProviderException(f"Incorrect e-mail address ({email_address}).")


class UseCache:
    """Used if Fasta Complementation Option is \"NCBI\"
       Sequences downloaded from NCBI are stored and reused
       reused between uses of this program."""
    pass


class FastaProvider(abc.ABC):
    """Pangraph building from maf file parameter. Ignored when -not_dag parameter is set.
       Maf file usually contains not full sequences but only parts of them, aligned to each other.
       To build an exact pangraph the full sequences must be retrieved from: ncbi or local file system."""

    @abc.abstractmethod
    def get_base(self, sequence_id: SequenceID, i: int):
        pass
