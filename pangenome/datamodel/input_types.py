import csv
from io import StringIO
from typing import Optional, Dict, Any
from pathlib import Path

from .Sequence import SequenceID


class InputError(Exception):
    pass


class Maf:
    """todo"""

    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        self.filename = filename
        self._raise_exception_if_incorrect(filecontent)
        self.filecontent = filecontent

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if False:
            raise InputError("Incorrect maf because...")

    @staticmethod
    def get_parameter_description():
        return 'Accepted formats: .maf, .po.'


class MetadataCSV:
    """todo See... examples\\Ebola\\ebola_metadata.csv"""

    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        filecontent = filecontent.read()
        self.filename: Path = filename
        self.metadata: Dict[SequenceID, Dict[str, Any]] = self.csv_to_dict(filecontent)

    @classmethod
    def csv_to_dict(cls, filecontent: str):
        rd = csv.DictReader(StringIO(filecontent), delimiter=',')
        cls._raise_exception_if_incorrect(rd)

        d = {}
        for row in rd:
            seqid = SequenceID(row['seqid'], skip_part_before_dot=False)
            d[seqid] = dict(row)
            del d[seqid]['seqid']
        return d

    @staticmethod
    def _raise_exception_if_incorrect(rd: csv.DictReader) -> None:
        headers = rd.fieldnames

        if 'seqid' not in headers:
            raise InputError('No seqid column in metadata csv.')

    def get_all_sequences_ids(self):
        return [*self.metadata.keys()]

    def get_sequence_metadata(self, seq_id):
        return self.metadata[seq_id]

class Po:
    """todo"""

    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        self.filename = filename
        self._raise_exception_if_incorrect(filecontent)
        self.filecontent = filecontent

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if False:
            raise InputError("Incorrect maf because...")


class MissingSymbol:
    """todo"""

    def __init__(self, symbol: str) -> 'MissingSymbol':
        if len(symbol) != 1:
            raise InputError('Missing symbol must be a single character.')
