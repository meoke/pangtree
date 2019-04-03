from pathlib import Path

from datamodel.fasta_providers.FastaProvider import FastaProvider


class FromFile(FastaProvider):
    def __init__(self, fastas_file: Path):
        pass