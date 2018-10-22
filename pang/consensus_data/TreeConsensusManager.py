from consensus_data.ConsensusNode import ConsensusNode
from graph.PathManager import PathManager
from .ConsensusTree import ConsensusTree
from .Errors import NoConsensus


class TreeConsensusManager(PathManager):
    def __init__(self, max_nodes_count):
        PathManager.__init__(self, max_nodes_count=max_nodes_count)
        self.consensus_tree = ConsensusTree()

    def get_root_node(self):
        return self.consensus_tree.get_root_node()

    def get_consensus(self, consensus_id):
        if consensus_id in self.path_names_to_array_id.values():
            return self.paths[consensus_id]
        raise NoConsensus

    def add_node(self, node: ConsensusNode, remapped_best_path):
        node_id = self.add_path(f"Consensus", remapped_best_path)
        node.consensus_id = node_id
        self.consensus_tree.add_node(node)
        return node_id

    def get_nodes(self):
        return self.consensus_tree.nodes

    def get_node(self, node_id):
        return self.consensus_tree.get_node(node_id)