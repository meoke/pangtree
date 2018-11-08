class ConsensusTree(object):
    def __init__(self):
        self.nodes = []

    def get_root_node(self):
        return self.nodes[0]

    def add_node(self, node):
        self.nodes.append(node)

    def get_node(self, node_id):
        return self.nodes[node_id]

    def remove_consensus_node(self, c_id):
        node_id_to_remove = [i for i, n in enumerate(self.nodes) if n.consensus_id == c_id][0]
        del self.nodes[node_id_to_remove]
