from pathlib import Path
from consensus.algorithm.FindCutoff import FindMaxCutoff, FindNodeCutoff
from metadata.MultialignmentMetadata import MultialignmentMetadata


class AlgorithmParams:
    outputdir: Path = None
    genomes_info: MultialignmentMetadata = None
    max_node_strategy: FindMaxCutoff = None
    node_cutoff_strategy: FindNodeCutoff = None
    stop: float = None
    re_consensus: bool = None
