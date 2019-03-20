from typing import NewType

from Bio import Entrez

from fasta_providers.FastaProvider import FastaProvider
from tools import loggingtools

EntrezSequenceID = NewType("EntrezSequenceID", str)

detailed_logger = loggingtools.get_logger("details")

class FromEntrezFastaProvider(FastaProvider):
    def __init__(self,):
        super().__init__()
        Entrez.email = "pedziadkiewicz@gmail.com"

    def get_source(self, sequenceID: EntrezSequenceID, start: int = None, end: int = None) -> str:
        detailed_logger.info(f"Downloading from entrez sequence {sequenceID}...")
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


