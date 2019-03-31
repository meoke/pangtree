import zipfile
from io import StringIO
from pathlib import Path
from typing import Dict

from Bio import SeqIO

from fasta_providers.FastaProvider import FastaProvider
from data.custom_types import SequenceID, Sequence


class FromFileFastaProvider(FastaProvider):
    def __init__(self, fastas_file: Path):
        super().__init__()
        self.fastas: Dict[SequenceID, Sequence] = self._read_fastas(fastas_file)

    def get_source(self, sequence_id: str, start: int = None, end: int = None):
        try:
            return self.fastas[sequence_id][start: end]
        except Exception:
            raise Exception(f"Cannot provide sequence with id: {sequence_id}. "
                            f"It was not included in provided fasta files.")

    def _read_fastas(self, fastas_file: Path) -> Dict[SequenceID, Sequence]:
        fastas_file_suffix = fastas_file.suffix
        if fastas_file_suffix == '':
            raise Exception("Cannot recognize fasta file format. Append extension to file name. "
                            "Available file formats: zip, fasta.")
        else:
            fastas_file_extensions = fastas_file_suffix[1:]
            if fastas_file_extensions == "zip":
                return self._read_zip(fastas_file)
            elif fastas_file_extensions == "fasta":
                return self._read_fasta(fastas_file)
            else:
                raise Exception("Unknown fasta file format. "
                                "Available file formats: zip, fasta, json.")

    def _read_zip(self, zip_path: Path):
        if not zipfile.is_zipfile(zip_path):
            raise Exception("Incorrect zip fasta source.")
        fastas_dict = {}
        with zipfile.ZipFile(zip_path) as myzip:
            zipinfolist = myzip.infolist()
            for zipinfo in zipinfolist:
                with myzip.open(zipinfo) as myfile:
                    fasta = myfile.read().decode("ASCII")
                    for record in SeqIO.parse(StringIO(fasta), "fasta"):
                        self._add_record_to_dict(record, fastas_dict)
        if len(fastas_dict) == 0:
            raise Exception("No sequences in zip provided as fasta source or incorrect fastas in zip.")
        return fastas_dict

    def _read_fasta(self, fasta_path: Path) -> Dict[SequenceID, Sequence]:
        fastas_dict = {}
        with open(fasta_path) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                self._add_record_to_dict(record, fastas_dict)
        if len(fastas_dict) == 0:
            raise Exception("No sequences in fasta provided as fasta source or incorrect fasta.")
        return fastas_dict

    def _add_record_to_dict(self, record, fastas_dict):
        if len(record.seq) == 0:
            raise Exception("Empty sequence in fasta source file. "
                            "Provide the sequence or remove its identifier.")
        if record.id in fastas_dict.keys():
            raise Exception("Incorrect fasta provided as fasta source. Sequences ids are not unique.")
        fastas_dict[FromFileFastaProvider.extract_seq_id(record.id)] = record.seq

    @staticmethod
    def extract_seq_id(fasta_id: str) -> SequenceID:
        splitted = fasta_id.split('.')
        if len(splitted) > 1:
            return ".".join(splitted[1:])
        elif len(splitted) == 1:
            return splitted[0]
