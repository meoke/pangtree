from pathlib import Path

from fasta_providers.FastaProvider import FastaProvider
from tools.pathtools import get_child_file_path

class FromZIPSystemProvider(FastaProvider):
    def __init__(self, fastas_dictionary: Path):
        super().__init__()
        self.fastas_dictionary = fastas_dictionary

    def get_source(self, sequenceID: str, start: int = None, end: int = None):
        fasta_path = get_child_file_path(self.fastas_dictionary, f"{sequenceID}.fasta")
        with open(fasta_path) as fasta:
            fasta_content = self.get_raw_sequence_from_fasta(fasta)
            return fasta_content[start:end]