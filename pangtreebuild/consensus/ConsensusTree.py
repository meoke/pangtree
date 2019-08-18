from typing import List, NewType, Dict, Optional

from newick import Node, dumps
from Bio import Phylo
from io import StringIO


from pangtreebuild.consensus.input_types import P
from pangtreebuild.datamodel.Sequence import SequenceID, SequencePath

ConsensusNodeID = NewType('ConsensusNodeID', int)


class ConsensusesTreeException(Exception):
    pass


class CompatibilityToPath:
    """Expresses how similar is one poagraph path to another."""

    def __init__(self, compatibility: float, p: P = P(1)):
        self.value: float = compatibility**p.value
        self.p: float = p.value

    def __eq__(self, other):
        return self.value == other.value

    def __lt__(self, other):
        return self.value < other.value

    def __le__(self, other):
        return self.value <= other.value

    def __gt__(self, other):
        return self.value > other.value

    def __ge__(self, other):
        return self.value >= other.value

    def __sub__(self, other):
        return CompatibilityToPath(self.value - other.value, P(1))

    def __str__(self):
        return f"""{self.value}"""

    def __repr__(self):
        return f"""{self.value}"""

    def root_value(self):
        return CompatibilityToPath(self.value**(1/self.p))


class ConsensusNode(object):
    def __init__(self,
                 consensus_id: ConsensusNodeID,
                 parent_node_id: Optional[ConsensusNodeID] = None,
                 children_nodes_ids: Optional[List[ConsensusNodeID]] = None,
                 sequences_ids: Optional[List[SequenceID]] = None,
                 mincomp: Optional[CompatibilityToPath] = None,
                 compatibilities_to_all: Optional[Dict[SequenceID, CompatibilityToPath]] = None,
                 consensus_path: Optional[SequencePath] = None):
        self.parent_node_id: ConsensusNodeID = parent_node_id
        self.children_nodes_ids: List[ConsensusNodeID] = children_nodes_ids if children_nodes_ids else []
        self.sequences_ids: List[SequenceID] = sequences_ids if sequences_ids else []
        self.consensus_id: ConsensusNodeID = consensus_id
        self.mincomp: CompatibilityToPath = mincomp if mincomp else CompatibilityToPath(0)
        self.compatibilities_to_all: Dict[SequenceID, CompatibilityToPath] = compatibilities_to_all if compatibilities_to_all else {}
        self.consensus_path: SequencePath = consensus_path

    def __str__(self):
        return f"ID: {self.consensus_id}, "\
            f"parentID: {self.parent_node_id}, " \
            f"children: {self.children_nodes_ids}, " \
            f"mincomp: {self.mincomp.value}, " \
            f"path length: {len(self.consensus_path)}, "\
            f"sequences ids: {self.sequences_ids}."

    def __eq__(self, other):
        return self.parent_node_id == other.parent_node_id and \
            self.consensus_id == other.consensus_id and \
            self.sequences_ids == other.sequences_ids and \
            abs(self.mincomp.value - other.mincomp.value) < 0.001 and \
            self.compatibilities_to_all == other.compatibilities_to_all and \
            self.consensus_path == other.consensus_path


class ConsensusTree:
    def __init__(self):
        self.nodes: List[ConsensusNode] = []

    def get_node(self, node_id: ConsensusNodeID) -> ConsensusNode:
        for node in self.nodes:
            if node.consensus_id == node_id:
                return node
        raise ConsensusesTreeException("No consensus with given ID.")

    def get_max_node_id(self):
        if len(self.nodes) == 0:
            return -1
        return max([node.consensus_id for node in self.nodes])

    def as_newick(self, seq_id_to_seq_name: Dict[SequenceID, str] = None, expand_leaves=False):
        return self._convert_to_newick(seq_id_to_seq_name, expand_leaves)

    def as_phyloxml(self, seq_id_to_seq_name: Dict[SequenceID, str] = None):
        return self._convert_to_phyloxml(seq_id_to_seq_name)

    def _convert_to_newick(self, seq_id_to_metadata: Dict[SequenceID, str] = None, expand_newick = False) -> str:
        def newick_nhx(newick_tree):

            label = newick_tree.name or ''
            if newick_tree._length:
                for cn in sorted_nodes:
                    if str(cn.consensus_id) == newick_tree.name:
                        if seq_id_to_metadata:
                            if len(cn.sequences_ids) == 1:
                                name = seq_id_to_metadata[cn.sequences_ids[0]]["name"]
                                group = seq_id_to_metadata[cn.sequences_ids[0]]["group"]
                                seqid = cn.sequences_ids[0]
                                metadata = f"[&&NHX:name={name}:group={group}:seqid={seqid}:mincomp={cn.mincomp}]"
                            elif len(cn.sequences_ids) == 0:
                                name = f"EmptyConsensus {cn.consensus_id}"
                                metadata = f"[&&NHX:name={name}:mincomp={cn.mincomp}]"
                            else:
                                name = f"Consensus {cn.consensus_id}"
                                metadata = f"[&&NHX:name={name}:mincomp={cn.mincomp}]"
                        else:
                            if len(cn.sequences_ids) == 1:
                                name = cn.sequences_ids[0]
                            elif len(cn.sequences_ids) == 0:
                                name = f"EmptyConsensus {cn.consensus_id}"
                            else:
                                name = f"Consensus {cn.consensus_id}"
                            mincomp = cn.mincomp
                            metadata = f"[&&NHX:name={name}:mincomp={mincomp}]"
                try:
                    label += ':' + newick_tree._length + metadata
                except:
                    print("metadata")
            descendants = ','.join([newick_nhx(n) for n in newick_tree.descendants])
            if descendants:
                descendants = '(' + descendants + ')'
            return descendants + label

        if not self.nodes:
            return None

        sorted_nodes = sorted(self.nodes, key=lambda x: x.consensus_id)

        if expand_newick:
            new_leaves_count = 0
            for node in self.nodes:
                if len(node.children_nodes_ids) == 0 and len(node.sequences_ids) > 1:
                    for seq_id in node.sequences_ids:
                        consensus_node_id = len(self.nodes) + new_leaves_count
                        node.children_nodes_ids.append(consensus_node_id)
                        sorted_nodes.append(ConsensusNode(consensus_id=consensus_node_id,
                                                           parent_node_id=node.consensus_id,
                                                           children_nodes_ids=[],
                                                           sequences_ids=[seq_id],
                                                           mincomp=CompatibilityToPath(1.0)
                                                           ))
                        new_leaves_count += 1

        nodes_to_process = [(None, sorted_nodes[0])]
        newick_tree = None
        while nodes_to_process:
            n = nodes_to_process.pop()
            node_parent_label = n[0]
            node = n[1]

            label = str(node.consensus_id)
            if node.parent_node_id is None:
                length = "1"
            else:
                parent_minComp = sorted_nodes[node.parent_node_id].mincomp.root_value().value
                length = str((1 - parent_minComp) - (1-node.mincomp.root_value().value))

            newick_node = Node(name=label, length=length)

            if newick_tree is None:
                newick_tree = newick_node
            else:
                parent_node = newick_tree.get_node(node_parent_label)
                parent_node.add_descendant(newick_node)

            for child in node.children_nodes_ids:
                nodes_to_process.append((label, sorted_nodes[child]))
        return "(" + newick_nhx(newick_tree) + ")"
