from pathlib import Path
from consensus.ConsensusesTree import ConsensusesTree
from datamodel.Pangraph import Pangraph
from tools import loggingtools

global_logger = loggingtools.get_global_logger()


class SimplePOAConsensusGenerator:
    def __init__(self, hbmin: float, blosum_path: Path):
        self.hbmin: float = hbmin
        self.blosum_path: Path = blosum_path

    def get_consensuses_tree(self,
                             pangraph: Pangraph,
                             output_dir: Path
                             ) -> ConsensusesTree:
        global_logger.info("Simple consensus generation started.")
        raise NotImplementedError()
