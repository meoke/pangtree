from typing import List, Dict, NewType
from pangraph.custom_types import SequenceID, NodeID

ConsensusNodeID = NewType('ConsensusNodeID', int)
Compatibility = NewType('Compatibility', float)


class ConsensusNode(object):
    def __init__(self,
                 parent_node_id: ConsensusNodeID = None,
                 children_nodes_ids: List[ConsensusNodeID] = None,
                 sequences_ids: List[SequenceID] = None,
                 consensus_id: ConsensusNodeID = None,
                 mincomp: Compatibility = None,
                 compatibilities_to_all: Dict[SequenceID, Compatibility] = None,
                 consensus_path: List[NodeID] = None):
        self.parent_node_id: ConsensusNodeID = parent_node_id
        self.children_nodes_ids: List[ConsensusNodeID] = children_nodes_ids if children_nodes_ids else []
        self.sequences_ids: List[SequenceID] = sequences_ids if sequences_ids else []
        self.consensus_id: ConsensusNodeID = consensus_id
        self.mincomp: Compatibility = mincomp
        self.compatibilities_to_all: Dict[SequenceID, Compatibility] = compatibilities_to_all
        self.consensus_path: List[NodeID] = consensus_path

    def __str__(self):
        return f"ID: {self.consensus_id}, "\
            f"parentID: {self.parent_node_id}, " \
            f"children: {self.children_nodes_ids}, " \
            f"mincomp: {self.mincomp}, " \
            f"path length: {len(self.consensus_path)}, "\
            f"sequences ids: {self.sequences_ids}."

    def __eq__(self, other):
        return self.parent_node_id == other.parent_node_id and \
            self.consensus_id == other.consensus_id and \
            self.sequences_ids == other.sequences_ids and \
            abs(self.mincomp - other.mincomp) < 0.001 and \
            self.compatibilities_to_all == other.compatibilities_to_all and \
            self.consensus_path == other.consensus_path
