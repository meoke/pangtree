from pathlib import Path
from consensus.ConsensusesTree import ConsensusesTree
from consensus.ConsensusNode import Compatibility
from pangraph.Pangraph import Pangraph
from metadata.MultialignmentMetadata import MultialignmentMetadata


class SimplePOAConsensusGenerator:
    def __init__(self, hbmin: Compatibility, blosum_path: Path):
        self.hbmin: Compatibility = hbmin
        self.blosum_path: Path = blosum_path

    def get_consensuses_tree(self,
                             pangraph: Pangraph,
                             genomes_info: MultialignmentMetadata,
                             output_dir: Path
                             ) -> ConsensusesTree:
        raise NotImplementedError()
