from pathlib import Path
from consensus.FindCutoff import FindMaxCutoff, FindNodeCutoff
from metadata.MultialignmentMetadata import MultialignmentMetadata


class AlgorithmParams:
    outputdir: Path = None
    cutoffs_log_path: Path = None
    genomes_info: MultialignmentMetadata = None
    max_node_strategy: FindMaxCutoff = None
    node_cutoff_strategy: FindNodeCutoff = None
    stop: float = None
    re_consensus: bool = None
