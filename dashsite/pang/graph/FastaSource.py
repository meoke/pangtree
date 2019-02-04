import abc
from Bio import Entrez


class FastaSource(abc.ABC):
    @abc.abstractmethod
    def get_source(self, id: str, start: int = None, end: int = None):
        pass


class EntrezFastaSource(FastaSource):
    def __init__(self,):
        super().__init__()
        Entrez.email = "pedziadkiewicz@gmail.com"

    def get_source(self, id: str, start: int = None, end: int = None) -> str:
        ncbi_id = id.split('.')[1]
        if 'v1' in ncbi_id:
            ncbi_id = ncbi_id.replace('v1', '.1')
        if 'v2' in ncbi_id:
            ncbi_id = ncbi_id.replace('v2', '.2')
        if 'v3' in ncbi_id:
            ncbi_id = ncbi_id.replace('v3', '.3')
        try:
            if start is not None and end is not None:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text", seq_start=start, seq_stop=end)

            else:
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
            fasta_first_line = handle.readline()
            fasta_content = handle.read().replace('\n', '')
            return fasta_content
        except Exception as ex:
            print(ncbi_id, start, end)
            raise ex


