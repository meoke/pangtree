import csv
from io import StringIO
from typing import Optional, Dict, Any, List
from pathlib import Path

from .Sequence import SequenceID


class InputError(Exception):
    pass


class Maf:
    """Multialignment content as stream. Accepted formats: maf and po."""

    def __init__(self, file_content: StringIO, filename: Optional[Path]):
        self.filename = filename
        self._raise_exception_if_incorrect(file_content)
        self.filecontent = file_content

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if False:
            raise InputError("Incorrect maf because...")


class MetadataCSV:
    """Metadata provides additional information about sequences used in Poagraph.
    See example metadata file: data\\Ebola\\ebola_metadata.csv"""

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
            if seqid in d:
                raise InputError("Not unique values seqid column in metadata file. Make them unique.")
            d[seqid] = dict(row)
            if None in d[seqid].keys():
                raise InputError("CSV metadata error. "
                                 "Different fields number in line 0 than in header line.")
            del d[seqid]['seqid']
        return d

    @staticmethod
    def _raise_exception_if_incorrect(rd: csv.DictReader) -> None:
        headers = rd.fieldnames
        if not headers:
            raise InputError('Empty csv file.')

        if 'seqid' not in headers:
            raise InputError('No \'seqid\' column in metadata csv.')

        if headers.count('seqid') > 1:
            raise InputError('Only one \'seqid\' column in metadata csv is allowed.')

    def get_all_sequences_ids(self) -> List[SequenceID]:
        return [*self.metadata.keys()]

    def get_sequence_metadata(self, seq_id) -> Dict[str, Any]:
        if seq_id in self.metadata:
            return self.metadata[seq_id]
        return {}

    def get_metadata_keys(self) -> List[str]:
        if self.metadata.values():
            example_metadata = [*self.metadata.values()][0]
            if example_metadata:
                return [*example_metadata.keys()]
        return []


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
    """A symbol used to complement omitted symbols in genome sequence. "?" is default."""

    def __init__(self, symbol: Optional[str] = '?'):
        if len(symbol) != 1:
            raise InputError('Missing symbol must be a single character.')
        self.value = symbol
