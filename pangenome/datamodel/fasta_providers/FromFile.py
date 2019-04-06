import zipfile
from pathlib import Path
from typing import Dict
from io import StringIO
from Bio import SeqIO

from datamodel.Node import Base
from datamodel.Sequence import SequenceID
from datamodel.fasta_providers.FastaProvider import FastaProvider, FastaProviderException


class FromFile(FastaProvider):
    def __init__(self, fastas_file: Path):
        self.sequences: Dict[SequenceID, str] = self._read_fastas(fastas_file)

    def get_base(self, sequence_id: SequenceID, i: int) -> Base:
        if sequence_id not in self.sequences.keys():
            raise FastaProviderException(f"Wrong sequence id: {sequence_id}. ")
        if i > len(self.sequences[sequence_id]):
            raise FastaProviderException(f"Index {i} is to large for sequence {SequenceID}.")
        return Base(self.sequences[sequence_id][i])

    def _read_fastas(self, fastas_file: Path) -> Dict[SequenceID, str]:
        fastas_file_suffix = fastas_file.suffix
        if fastas_file_suffix == '':
            raise FastaProviderException("No file extension in fasta source. Append extension to file name."
                                         "Available file formats: zip, fasta, fna, faa.")
        else:
            fastas_file_extensions = fastas_file_suffix[1:]
            if fastas_file_extensions == "zip":
                return self._read_zip(fastas_file)
            elif fastas_file_extensions in ["fna", "fasta", "faa"]:
                return self._read_fasta(fastas_file)
            else:
                raise FastaProviderException("Unknown fasta file format. "
                                             "Available file formats: zip, fasta, json.")

    def _read_zip(self, zip_path: Path) -> Dict[SequenceID, str]:
        if not zipfile.is_zipfile(zip_path):
            raise FastaProviderException("Incorrect zip fasta source.")
        fastas_dict = {}
        with zipfile.ZipFile(zip_path) as myzip:
            zipinfolist = myzip.infolist()
            for zipinfo in zipinfolist:
                with myzip.open(zipinfo) as myfile:
                    fasta = myfile.read().decode("ASCII")
                    for record in SeqIO.parse(StringIO(fasta), "fasta"):
                        self._add_record_to_dict(record, fastas_dict)
        if len(fastas_dict) == 0:
            raise FastaProviderException("No sequences in zip provided as fasta source or incorrect fastas in zip.")
        return fastas_dict

    def _read_fasta(self, fasta_path: Path) -> Dict[SequenceID, str]:
        fastas_dict = {}
        with open(fasta_path) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                self._add_record_to_dict(record, fastas_dict)
        if len(fastas_dict) == 0:
            raise FastaProviderException("No sequences in fasta provided as fasta source or incorrect fasta.")
        return fastas_dict

    def _add_record_to_dict(self, record, fastas_dict):
        if len(record.seq) == 0:
            raise FastaProviderException("Empty sequence in fasta source file. "
                                         "Provide the sequence or remove its identifier.")
        if SequenceID(record.id) in fastas_dict.keys():
            raise FastaProviderException("Incorrect fasta provided as fasta source. Sequences ids are not unique.")
        fastas_dict[SequenceID(record.id)] = record.seq
