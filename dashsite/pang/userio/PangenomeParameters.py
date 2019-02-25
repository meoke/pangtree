from typing import Optional

from io import StringIO
from enum import Enum
from pathlib import Path


class ConsensusAlgorithm(Enum):
    NO = 0
    SIMPLE = 1
    TREE = 2


class FastaComplementationOption(Enum):
    NO = 0
    NCBI = 1
    LOCAL = 2


class MaxCutoffOption(Enum):
    MAX1 = 1
    MAX2 = 2


class NodeCutoffOption(Enum):
    NODE1 = 1
    NODE2 = 2
    NODE3 = 3
    NODE4 = 4


class PangenomeParameters:
    def __init__(self,
                 multialignment_file_content: StringIO,
                 multialignment_file_path: Path,
                 metadata_file_content: Optional[str],
                 metadata_file_path: Path,
                 output_path: Path,
                 generate_fasta: bool,
                 consensus_type: ConsensusAlgorithm,
                 hbmin: float,
                 range: Optional[list],
                 multiplier: Optional[float],
                 stop: Optional[float],
                 re_consensus: Optional[bool],
                 not_dag: bool,
                 fasta_complementation_option: Optional[FastaComplementationOption],
                 local_fasta_dirpath: Optional[Path],
                 max_cutoff_option: Optional[MaxCutoffOption],
                 node_cutoff_option: Optional[NodeCutoffOption]
                 ):
        self.multialignment_file_content = multialignment_file_content
        self.multialignment_file_path = multialignment_file_path
        self.metadata_file_content = metadata_file_content
        self.metadata_file_path = metadata_file_path
        self.output_path = output_path

        self.not_dag = not_dag
        self.fasta_complementation_option = fasta_complementation_option
        self.local_fasta_dirpath = local_fasta_dirpath

        self.generate_fasta = generate_fasta

        self.consensus_type = consensus_type

        self.hbmin = hbmin

        self.stop = stop
        self.max_cutoff_option = max_cutoff_option
        self.range = range
        self.node_cutoff_option = node_cutoff_option
        self.multiplier = multiplier
        self.re_consensus = re_consensus

        self._validate()

    def _validate(self):
        if self.multialignment_file_content is None:
            raise Exception("Unspecified multialignment file.")

        if self.output_path is None:
            raise Exception("Unspecified output path.")

        if self.fasta_complementation_option is FastaComplementationOption.LOCAL and self.local_fasta_dirpath is None:
            raise Exception("Unspecified path to direction with fasta files, "
                            "while FastaComplementationOption.LocalFasta was chosen. "
                            "Use FastaComplementationOption.No or FastaComplementationOption.NCBI instead.")

        if self.consensus_type is None:
            raise Exception("Unspecified consensus algorithm. "
                            "Use ConsensusAlgorithm.No if no consensus generation should be generated.")

        if self.consensus_type is ConsensusAlgorithm.SIMPLE and not self.hbmin:
            raise Exception("For SIMPLE consensus algorithm HBMIN must be specified.")

        if self.consensus_type is ConsensusAlgorithm.TREE:
            if self.stop is None:
                raise Exception("For TREE consensus algorithm STOP must be specified.")
            if self.max_cutoff_option is None:
                raise Exception("For TREE consensus algorithm MAX CUTOFF OPTION must be specified.")
            if self.node_cutoff_option is None:
                raise Exception("For TREE consensus algorithm NODE CUTOFF OPTION must be specified.")
            if self.re_consensus is None:
                raise Exception("For TREE consensus algorithm RE CONSENSUS must be specified.")

        if self.max_cutoff_option is MaxCutoffOption.MAX1:
            if self.range is None:
                raise Exception("For MAX1 max cutoff option CUTOFF SEARCH RANGE must be specified.")
            if len(self.range) != 2:
                raise Exception("CUTOFF SEARCH RANGE must have length 2.")
            if self.range[1] < self.range[0]:
                raise Exception("CUTOFF SEARCH RANGE first value must be smaller or equal to second value.")
            if self.range[0] < 0 or self.range[0] > 1 or self.range[1] < 0 or self.range[1] > 1:
                raise Exception("CUTOFF SEARCH RANGE values must be in the range of [0,1].")

        if self.node_cutoff_option in [NodeCutoffOption.NODE1, NodeCutoffOption.NODE2]:
            if self.multiplier is None:
                raise Exception("For NODE1 and NODE2 node cutoff option MULTIPLIER must be specified.")

        if self.stop < 0 or self.stop > 1:
            raise Exception("STOP value must be in the range of [0,1].")





