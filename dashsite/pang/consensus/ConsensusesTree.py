from consensus.errors import TreeConsensusGenerationException


class ConsensusesTree:
    def __init__(self):
        self.nodes = []

    def get_node(self, node_id):
        for node in self.nodes:
            if node.consensus_id == node_id:
                return node
        raise TreeConsensusGenerationException("No consensus with given ID.")