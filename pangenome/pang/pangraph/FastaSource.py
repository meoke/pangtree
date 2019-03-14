import abc
from pathlib import Path
from typing import NewType

from Bio import Entrez

from pangraph.custom_types import SequenceID
from tools.pathtools import get_child_file_path

EntrezSequenceID = NewType("EntrezSequenceID", str)

class FastaSource(abc.ABC):
    @abc.abstractmethod
    def get_source(self, sequenceID: str, start: int = None, end: int = None):
        pass

    def get_raw_sequence_from_fasta(self, fasta_handle):
        _ = fasta_handle.readline()
        return fasta_handle.read().replace('\n', '')


class FastaFileSystemSource(FastaSource):
    def __init__(self, fastas_dictionary: Path):
        super().__init__()
        self.fastas_dictionary = fastas_dictionary

    def get_source(self, sequenceID: str, start: int = None, end: int = None):
        fasta_path = get_child_file_path(self.fastas_dictionary, f"{sequenceID}.fasta")
        with open(fasta_path) as fasta:
            fasta_content = self.get_raw_sequence_from_fasta(fasta)
            return fasta_content[start:end]


class EntrezFastaSource(FastaSource):
    def __init__(self,):
        super().__init__()
        Entrez.email = "pedziadkiewicz@gmail.com"

    def get_source(self, sequenceID: EntrezSequenceID, start: int = None, end: int = None) -> str:
        try:
            if start is not None and end is not None:
                handle = Entrez.efetch(db="nucleotide",
                                       id=sequenceID,
                                       rettype="fasta",
                                       retmode="text",
                                       seq_start=start,
                                       seq_stop=end)
            else:
                handle = Entrez.efetch(db="nucleotide", id=sequenceID, rettype="fasta", retmode="text")
            fasta_content = self.get_raw_sequence_from_fasta(handle)
            return fasta_content
        except Exception as ex:
            raise Exception(f"Cannot download from Entrez sequence of ID: {sequenceID}") from ex


