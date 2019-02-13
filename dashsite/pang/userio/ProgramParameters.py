from enum import Enum
from pathlib import Path


class ConsensusAlgorithm(Enum):
    No = 0
    Simple = 1
    Tree = 2


class FastaComplementationOption(Enum):
    No = 0
    NCBI = 1
    LocalFasta = 2

class ProgramParameters:
    def __init__(self):
        self.multialignment_file_path: Path = None
        self.metadata_file_path: Path = None
        self.output_path: Path = None
        self.generate_fasta: bool = None
        self.consensus_type: ConsensusAlgorithm = None
        self.hbmin: float = None
        self.r: list = None
        self.multiplier: float = None
        self.stop: float = None
        self.re_consensus: bool = None
        self.anti_granular: bool = None
        self.not_dag: bool = None
        self.fasta_complementation: FastaComplementationOption = None
        self.local_fasta_dirpath: Path = None
