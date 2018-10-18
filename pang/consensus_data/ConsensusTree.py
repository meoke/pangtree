from .ConsensusNode import ConsensusNode


class ConsensusTree(object):
    def __init__(self):
        self.nodes = [ConsensusNode(consensus_id=0)]

    def get_root_node(self):
        return self.nodes[0]

    def add_node(self, node):
        self.nodes.append(node)