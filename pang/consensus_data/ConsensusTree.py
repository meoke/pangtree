from .ConsensusNode import ConsensusNode


class ConsensusTree(object):
    def __init__(self):
        self.nodes = [ConsensusNode(consensus_id=0)]

    def get_root_node(self):
        return self.nodes[0]