import abc
from Bio import Entrez


class FastaSource(abc.ABC):
    @abc.abstractmethod
    def get_source(self, id: str, start: int, end: int):
        pass


class EntrezFastaSource(FastaSource):
    def __init__(self, email):
        Entrez.email = email

    def get_source(self, id: str, start: int, end: int) -> str:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)
        return handle.read()
