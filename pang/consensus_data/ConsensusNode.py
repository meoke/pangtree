class ConsensusNode(object):
    def __init__(self, parent_node_id=None, children_nodes=None, sequences_names=None, consensus_id=-1, mincomp=None):
        self.children_nodes = children_nodes if children_nodes else []
        self.sequences_names = sequences_names if sequences_names else []
        self.consensus_id = consensus_id
        self.parent_node_id = parent_node_id
        self.compatibilities_to_all = None
        self.mincomp = mincomp

    def get_compatibilities_to_own_sources(self):
        return [self.compatibilities_to_all[seqname] for seqname in self.sequences_names]
