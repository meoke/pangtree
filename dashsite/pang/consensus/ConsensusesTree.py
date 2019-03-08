from typing import List
from consensus.ConsensusNode import ConsensusNodeID, ConsensusNode
from consensus.exceptions import ConsensusesTreeException


class ConsensusesTree:
    def __init__(self):
        self.nodes: List[ConsensusNode] = []

    def get_node(self, node_id: ConsensusNodeID) -> ConsensusNode:
        for node in self.nodes:
            if node.consensus_id == node_id:
                return node
        raise ConsensusesTreeException("No consensus with given ID.")
