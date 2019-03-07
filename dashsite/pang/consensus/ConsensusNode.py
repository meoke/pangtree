from typing import List, Dict
from graph2.custom_types import SequenceID, Sequence


class ConsensusNode(object):
    def __init__(self,
                 parent_node_id: int = None,
                 children_nodes_ids: List[int] = None,
                 sequences_ids: List[SequenceID] = None,
                 consensus_id: int = None,
                 mincomp=None,
                 compatibilities_to_all: Dict[SequenceID, float] = None,
                 consensus_path: Sequence = None):
        self.parent_node_id = parent_node_id
        self.children_nodes_ids = children_nodes_ids if children_nodes_ids else []
        self.sequences_ids = sequences_ids if sequences_ids else []
        self.consensus_id = consensus_id
        self.mincomp = mincomp
        self.compatibilities_to_all = compatibilities_to_all
        self.consensus_path = consensus_path

    # def get_compatibilities_to_own_sources(self):
    #     return [self.compatibilities_to_all[seq_id] for seq_id in self.sequences_ids]
