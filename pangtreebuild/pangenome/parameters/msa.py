"""Objects to store multiple sequence alignment files content: MAF, PO."""

import csv
from io import StringIO
from pathlib import Path
from typing import Optional, Any, Dict, List


class Maf:
    """MSA in MAF format. https://genome.ucsc.edu/FAQ/FAQformat.html#format5

    Args:
        file_name: MAF file name, no requirement to have full path.
        file_content: MAF content as stream.

    Raises:
        ValueError: If MAF is incorrect.

    Attributes:
        file_name (Path): File name.
        filecontent (str): MSA (MAF or PO) content content as stream."""

    def __init__(self, file_content: StringIO, file_name: Optional[Path]):
        self.filename = file_name
        self._raise_exception_if_incorrect(file_content)
        self.filecontent = file_content

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if not filecontent:
            raise ValueError("Incorrect MAF - empty file.")


class Po:
    """Multialignment in PO format.

    https://github.com/meoke/pangtree/blob/master/Documentation.md#po-file-format-specification

    Args:
        file_name: MAF file name, no requirement to have full path.
        file_content: MAF content as stream.

    Raises:
        ValueError: If MAF is incorrect.

    Attributes:
        file_name (Path): File name.
        file_content (str): MSA (MAF or PO) content content as stream."""

    def __init__(self, file_content: StringIO, file_name: Optional[Path]):
        self.filename = file_name
        self._raise_exception_if_incorrect(file_content)
        self.filecontent = file_content

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if not filecontent:
            raise ValueError("Incorrect PO - empty file.")


class SequenceID(object):
    """ID of genome seuquence being a part of poagraph.

    Args:
        sequence_id: value of the sequence ID as str
        skip_part_before_dot: TODO remove this hack so that
                            always part before dot is taken as sequence ID

    Attributes:
        value (str): sequence ID
    """

    value: str

    def __init__(self, sequence_id: str):
        if len(sequence_id) == 0:
            raise ValueError("Sequence ID cannot be empty.")
        if sequence_id == "eboVir3.KM034562v1":  # Ebola data hack
            self.value = "KM034562v1"
        else:
            self.value = sequence_id.split('.')[0]

    def __str__(self):
        return f"{self.value}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other: 'SequenceID'):
        return other and self.value.__eq__(other.value)

    def __hash__(self):
        return hash(self.value)

    def __le__(self, other: 'SequenceID'):
        return self.value.__le__(other.value)

    def __lt__(self, other: 'SequenceID'):
        return self.value.__lt__(other.value)

    def __ge__(self, other: 'SequenceID'):
        return self.value.__ge__(other.value)

    def __gt__(self, other: 'SequenceID'):
        return self.value.__gt__(other.value)


class MetadataCSV:
    """Metadata provides additional information about sequences.
    See example metadata file: data\\Ebola\\ebola_metadata.csv

    Args:
        filecontent: Content of the CSV file with metadata.
        filename: Name of the file.

    Attributes:
        metadata: Dictionary .
        filename: Name of the file.
    """

    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        filecontent = filecontent.read()
        self.filename: Path = filename
        self.metadata: Dict[SequenceID, Dict[str, Any]]
        self.metadata = self.csv_to_dict(filecontent)

    @classmethod
    def csv_to_dict(cls, filecontent: str):
        rd = csv.DictReader(StringIO(filecontent), delimiter=',')
        cls._raise_exception_if_incorrect(rd)

        d = {}
        for row in rd:
            seqid = SequenceID(row['seqid'])
            if seqid in d:
                raise ValueError("""Repeated values in seqid column in metadata file. Make them unique.""")
            d[seqid] = dict(row)
            if None in d[seqid].keys():
                raise ValueError("""CSV metadata error. Different number of columns in line 0 than in header line.""")
            del d[seqid]['seqid']
        return d

    @staticmethod
    def _raise_exception_if_incorrect(rd: csv.DictReader) -> None:
        headers = rd.fieldnames
        if not headers:
            raise ValueError('Empty csv file.')

        if 'seqid' not in headers:
            raise ValueError('No \'seqid\' column in metadata csv.')

        if headers.count('seqid') > 1:
            raise ValueError("""Only one \'seqid\' column in metadata csv is allowed.""")

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
