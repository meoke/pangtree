from pathlib import Path
from consensus.ConsensusesTree import ConsensusesTree
from consensus.ConsensusNode import Compatibility
from pangraph.Pangraph import Pangraph


class SimplePOAConsensusGenerator:
    def __init__(self, hbmin: Compatibility, blosum_path: Path):
        self.hbmin: Compatibility = hbmin
        self.blosum_path: Path = blosum_path

    def get_consensuses_tree(self,
                             pangraph: Pangraph,
                             output_dir: Path
                             ) -> ConsensusesTree:
        raise NotImplementedError()
