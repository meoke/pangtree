from graph.PathManager import PathManager
from .ConsensusTree import ConsensusTree
from .Errors import NoConsensus


class TreeConsensusManager(PathManager):
    def __init__(self):
        PathManager.__init__(self)
        self.consensus_tree = ConsensusTree()

    def get_root_node(self):
        return self.consensus_tree.get_root_node()

    def get_consensus(self, consensus_id):
        if consensus_id in self.path_names_to_array_id.keys():
            return self.paths[consensus_id]
        raise NoConsensus