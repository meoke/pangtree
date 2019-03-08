from typing import List, Dict, NewType
from pangraph.custom_types import SequenceID, Sequence

ConsensusNodeID = NewType('ConsensusNodeID', int)
Compatibility = NewType('Compatibility', float)

class ConsensusNode(object):
    def __init__(self,
                 parent_node_id: ConsensusNodeID = None,
                 children_nodes_ids: List[ConsensusNodeID] = None,
                 sequences_ids: List[SequenceID] = None,
                 consensus_id: ConsensusNodeID = None,
                 mincomp: Compatibility=None,
                 compatibilities_to_all: Dict[SequenceID, Compatibility] = None,
                 consensus_path: Sequence = None):
        self.parent_node_id: ConsensusNodeID = parent_node_id
        self.children_nodes_ids: List[ConsensusNodeID] = children_nodes_ids if children_nodes_ids else []
        self.sequences_ids: List[SequenceID] = sequences_ids if sequences_ids else []
        self.consensus_id: ConsensusNodeID = consensus_id
        self.mincomp: Compatibility = mincomp
        self.compatibilities_to_all: Dict[SequenceID, Compatibility] = compatibilities_to_all
        self.consensus_path: Sequence = consensus_path

