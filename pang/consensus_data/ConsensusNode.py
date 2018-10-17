class ConsensusNode(object):
    def __init__(self, children_nodes=None, sequences_ids=None, consensus_id=-1):
        self.children_nodes = children_nodes if children_nodes else []
        self.sequences_ids = sequences_ids if sequences_ids else []
        self.consensus_id = consensus_id