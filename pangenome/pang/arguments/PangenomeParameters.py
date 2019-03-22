import os
from typing import Optional

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
    def __init__(self, multialignment_file_content: str, multialignment_file_path: Path,
                 metadata_file_content: Optional[str], metadata_file_path: Optional[Path],
                 blosum_file_path: Optional[Path], output_path: Path, generate_fasta: bool,
                 consensus_type: ConsensusAlgorithm, hbmin: float, search_range: Optional[list],
                 multiplier: Optional[float], stop: Optional[float], re_consensus: Optional[bool], not_dag: bool,
                 fasta_complementation_option: Optional[FastaComplementationOption],
                 missing_nucleotide_symbol: Optional[str], local_fasta_dirpath: Optional[Path],
                 max_cutoff_option: Optional[MaxCutoffOption], node_cutoff_option: Optional[NodeCutoffOption],
                 verbose: bool, quiet: bool, email_address: str, cache: bool, p: float):
        self.multialignment_file_content = multialignment_file_content
        self.multialignment_file_path = multialignment_file_path
        self.metadata_file_content = metadata_file_content
        self.metadata_file_path = metadata_file_path
        self.blosum_file_path = blosum_file_path or self._get_default_blosum_path()
        self.output_path = output_path

        self.not_dag = not_dag
        self.fasta_complementation_option = fasta_complementation_option
        self.email_address = email_address
        self.cache = cache
        self.missing_nucleotide_symbol = missing_nucleotide_symbol or '?'
        self.local_fasta_dirpath = local_fasta_dirpath

        self.generate_fasta = generate_fasta

        self.consensus_type = consensus_type

        self.hbmin = hbmin

        self.stop = stop
        self.max_cutoff_option = max_cutoff_option
        self.search_range = search_range
        self.node_cutoff_option = node_cutoff_option
        self.multiplier = multiplier
        self.re_consensus = re_consensus
        self.verbose = verbose
        self.quiet = quiet
        self.p = p

        self._validate()

    def _validate(self) -> None:
        if self.multialignment_file_content is None:
            raise Exception("Unspecified multialignment file.")

        if self.output_path is None:
            raise Exception("Unspecified output path.")

        if self.fasta_complementation_option is FastaComplementationOption.NO and \
                self.missing_nucleotide_symbol is not None:
            try:
                self._blosum_contains_missing_nucl_symbol(blosum_path= self.blosum_file_path,
                                                          missing_nucleotide_symbol = self.missing_nucleotide_symbol)
            except Exception as e:
                raise Exception(f"The BLOSUM does not contain "
                                f"symbol specified for missing nucleotides ({self.missing_nucleotide_symbol}.") from e
        elif self.fasta_complementation_option is FastaComplementationOption.NO and \
            self.missing_nucleotide_symbol is None:
            try:
                self._blosum_contains_missing_nucl_symbol(missing_nucleotide_symbol = self.missing_nucleotide_symbol)
            except Exception as e:
                raise Exception("The BLOSUM does not contain default symbol for missing nucleotides (\'?\').") from e

        if self.fasta_complementation_option is FastaComplementationOption.LOCAL and self.local_fasta_dirpath is None:
            raise Exception("Unspecified path to direction with fasta files, "
                            "while FastaComplementationOption.LocalFasta was chosen. "
                            "Use FastaComplementationOption.No or FastaComplementationOption.NCBI instead.")

        if self.fasta_complementation_option is FastaComplementationOption.NCBI and self.email_address is None:
            raise Exception("Unspecified email address. "
                            "Email address is requiered if FastaComplementationOption.NCBI was chosen.")

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
            if self.search_range is None:
                raise Exception("For MAX1 max cutoff option CUTOFF SEARCH RANGE must be specified.")
            if len(self.search_range) != 2:
                raise Exception("CUTOFF SEARCH RANGE must have length 2.")
            if self.search_range[1] < self.search_range[0]:
                raise Exception("CUTOFF SEARCH RANGE first value must be smaller or equal to second value.")
            if self.search_range[0] < 0 \
                    or self.search_range[0] > 1 \
                    or self.search_range[1] < 0\
                    or self.search_range[1] > 1:
                raise Exception("CUTOFF SEARCH RANGE values must be in the range of [0,1].")

        if self.node_cutoff_option in [NodeCutoffOption.NODE1, NodeCutoffOption.NODE2]:
            if self.multiplier is None:
                raise Exception("For NODE1 and NODE2 node cutoff option MULTIPLIER must be specified.")

        if self.stop < 0 or self.stop > 1:
            raise Exception("STOP value must be in the range of [0,1].")

    def _blosum_contains_missing_nucl_symbol(self,
                                             blosum_path: Path,
                                             missing_nucleotide_symbol: str):
        with open(blosum_path) as b:
            blosum_lines = b.readlines()
        for blosum_line in blosum_lines:
            if len(blosum_line) > 0 and blosum_line[0] == " ":
                blosum_symbols_line = blosum_line
                break
        blosum_symbols = str.strip(blosum_symbols_line).split(" ")
        if missing_nucleotide_symbol in blosum_symbols:
            return True
        else:
            raise Exception("Error while a try of finding symbol for missing nucletides.")

    def _get_default_blosum_path(self):
        return Path(os.path.abspath(__file__)).joinpath('../../bin/blosum80.mat').resolve()

    def __str__(self):
        return f"""
        Pangenome parameters:
        multialignment_file_path: {self.multialignment_file_path}
        metadata_file_path: {self.metadata_file_path}
        blosum_file_path: {self.blosum_file_path}
        output_path: {self.output_path}
        not_dag: {self.not_dag}
        fasta_complementation_option: {self.fasta_complementation_option}
        missing_nucleotide_symbol: {self.missing_nucleotide_symbol}
        local_fasta_dirpath: {self.local_fasta_dirpath}
        generate_fasta: {self.generate_fasta}
        consensus_type: {self.consensus_type}
        hbmin: {self.hbmin}
        stop: {self.stop}
        max_cutoff_option: {self.max_cutoff_option}
        search_range: {self.search_range}
        node_cutoff_option: {self.node_cutoff_option}
        multiplier: {self.multiplier}
        re_consensus: {self.re_consensus},
        quiet: {self.quiet},
        verbose: {self.verbose},
        e-mail: {self.email_address},
        cache: {self.cache},
        p: {self.p}"""
