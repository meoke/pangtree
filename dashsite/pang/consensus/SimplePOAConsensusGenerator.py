from pathlib import Path
from collections import deque
from typing import List, Dict, Tuple

from consensus.ConsensusesTree import ConsensusesTree
from consensus.ConsensusNode import ConsensusNode, Compatibility
from consensus.FindCutoff import FindMaxCutoff, FindNodeCutoff
from consensus.exceptions import TreeConsensusGenerationException
from pangraph.Pangraph import Pangraph
from pangraph.custom_types import SequenceID, Sequence
from metadata.MultialignmentMetadata import MultialignmentMetadata
import consensus.top_consensus as top_consensus


class SimplePOAConsensusGenerator:
    def __init__(self,
                 hbmin: Compatibility):
        self.hbmin: Compatibility = hbmin

    def get_consensuses_tree(self,
                             pangraph: Pangraph,
                             genomes_info: MultialignmentMetadata,
                             output_dir: Path
                             ) -> ConsensusesTree:
        raise NotImplementedError()
