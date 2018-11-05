class ConsensusTree(object):
    def __init__(self):
        self.nodes = []

    def get_root_node(self):
        return self.nodes[0]

    def add_node(self, node):
        self.nodes.append(node)

    def get_node(self, node_id):
        return self.nodes[node_id]
