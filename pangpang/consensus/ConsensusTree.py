from typing import List, NewType, Dict

from consensus.input_types import P
from datamodel.Node import NodeID
from datamodel.Sequence import SequenceID, SequencePath

ConsensusNodeID = NewType('ConsensusNodeID', int)


class ConsensusesTreeException(Exception):
    pass


class CompatibilityToPath:
    """Expresses how similar is one poagraph path to another."""

    def __init__(self, compatibility: float, p: P):
        self.value: float = compatibility**p.value

    def __eq__(self, other):
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __le__(self, other):
        return self.value <= other.value

    def __gt__(self, other):
        return self.value > other.value

    def __ge__(self, other):
        return self.value >= other.value

    def __sub__(self, other):
        return CompatibilityToPath(self.value - other.value, P(1))

    def __str__(self):
        return f"""{self.value}"""

    def __repr__(self):
        return f"""{self.value}"""


class ConsensusNode(object):
    def __init__(self,
                 parent_node_id: ConsensusNodeID = None,
                 children_nodes_ids: List[ConsensusNodeID] = None,
                 sequences_ids: List[SequenceID] = None,
                 consensus_id: ConsensusNodeID = None,
                 mincomp: CompatibilityToPath = None,
                 compatibilities_to_all: Dict[SequenceID, CompatibilityToPath] = None,
                 consensus_path: SequencePath = None):
        self.parent_node_id: ConsensusNodeID = parent_node_id
        self.children_nodes_ids: List[ConsensusNodeID] = children_nodes_ids if children_nodes_ids else []
        self.sequences_ids: List[SequenceID] = sequences_ids if sequences_ids else []
        self.consensus_id: ConsensusNodeID = consensus_id
        self.mincomp: CompatibilityToPath = mincomp
        self.compatibilities_to_all: Dict[SequenceID, CompatibilityToPath] = compatibilities_to_all
        self.consensus_path: SequencePath = consensus_path

    def __str__(self):
        return f"ID: {self.consensus_id}, "\
            f"parentID: {self.parent_node_id}, " \
            f"children: {self.children_nodes_ids}, " \
            f"mincomp: {self.mincomp.value}, " \
            f"path length: {len(self.consensus_path)}, "\
            f"sequences ids: {self.sequences_ids}."

    def __eq__(self, other):
        return self.parent_node_id == other.parent_node_id and \
            self.consensus_id == other.consensus_id and \
            self.sequences_ids == other.sequences_ids and \
            abs(self.mincomp.value - other.mincomp.value) < 0.001 and \
            self.compatibilities_to_all == other.compatibilities_to_all and \
            self.consensus_path == other.consensus_path


class ConsensusTree:
    def __init__(self):
        self.nodes: List[ConsensusNode] = []

    def get_node(self, node_id: ConsensusNodeID) -> ConsensusNode:
        for node in self.nodes:
            if node.consensus_id == node_id:
                return node
        raise ConsensusesTreeException("No consensus with given ID.")

    def get_max_node_id(self):
        if len(self.nodes) == 0:
            return -1
        return max([node.consensus_id for node in self.nodes])
