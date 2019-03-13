import abc
from pathlib import Path

from Bio import Entrez

from pangraph.custom_types import SequenceID
from tools.pathtools import get_child_file_path


class FastaSource(abc.ABC):
    @abc.abstractmethod
    def get_source(self, sequenceID: SequenceID, start: int = None, end: int = None):
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

    def get_source(self, sequenceID: str, start: int = None, end: int = None) -> str:
        ncbi_id = sequenceID
        if 'v1' in ncbi_id:
            ncbi_id = ncbi_id.replace('v1', '.1')
        if 'v2' in ncbi_id:
            ncbi_id = ncbi_id.replace('v2', '.2')
        if 'v3' in ncbi_id:
            ncbi_id = ncbi_id.replace('v3', '.3')
        try:
            if start is not None and end is not None:
                handle = Entrez.efetch(db="nucleotide",
                                       id=ncbi_id,
                                       rettype="fasta",
                                       retmode="text",
                                       seq_start=start,
                                       seq_stop=end)
            else:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
            fasta_content = self.get_raw_sequence_from_fasta(handle)
            return fasta_content
        except Exception as ex:
            print(ncbi_id, start, end)
            raise ex


