from typing import NewType

from Bio import Entrez
from Bio import SeqIO

from fasta_providers.FastaProvider import FastaProvider
from tools import loggingtools, pathtools

EntrezSequenceID = NewType("EntrezSequenceID", str)

detailed_logger = loggingtools.get_logger("details")

class FromEntrezFastaProvider(FastaProvider):
    def __init__(self, email_address):
        super().__init__()
        Entrez.email = email_address

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


class FastaCache:
    def __init__(self, parent_dir):
        self.parent_dir = parent_dir
        self.cache_dir = pathtools.get_child_path(parent_dir, ".fastacache")

    def cache_dir_exists(self):
        return pathtools.dir_exists(self.cache_dir)

    def create_cahce_dir(self) -> None:
        pathtools.create_dir(self.cache_dir)

    def save_to_cache(self, seq_id, sequence)-> None:
        cache_filename = f"{seq_id}.fasta"
