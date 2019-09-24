from typing import List, NewType, Dict, Optional
import math
from newick import Node
from pangtreebuild.affinitytree.input_types import P
from pangtreebuild.datamodel.Sequence import SequenceID, SequencePath

AffinityNodeID = NewType('AffinityNodeID', int)


class AffinityTreeException(Exception):
    """Any exception connected with Affinity Tree"""

    pass


class Compatibility:
    """Asymetric similiarity measure of two poagraph paths.

    Attributes:
        compatibility (float): Raw compatibility value - count of common nodes devided by length of one of the paths.
        p (P): Parameter to control compatibility value interpretation. Compatibility is raised to the power of P.
    """

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
        if self.p != other.p:
            raise AffinityTreeException("Cannot subtract compatibilities: different P values.")
        return Compatibility(self.value - other.value, P(self.p))

    def __str__(self):
        return f"""{self.value}"""

    def __repr__(self):
        return f"""{self.value}"""

    def base_value(self):
        """Get compatibility value without transformation using P."""
        return Compatibility(self.value ** (1 / self.p))


class AffinityNode(object):
    """Affinity Tree node

    Attributes:
        id (AffinityNodeID): Affinity Node ID.
        parent (AffinityNodeID): ID of the parent node.
        children(List[AffinityNodeID]): IDs of the children nodes.
        sequences(List[SequenceID]): IDs of the sequences assigned to this node.
        mincomp(Compatibility): minimum from the compatibilities of this consensus to the assigned sequences
        compatibilities(Dict[SequenceID, Compatibility]): Dictionary of compatibilities to any sequences.
        consensus(SequencePath): path of the consensus defined as path in corresponding Poagraph.
    """

    def __init__(self,
                 id: AffinityNodeID,
                 parent: Optional[AffinityNodeID] = None,
                 children: Optional[List[AffinityNodeID]] = None,
                 sequences: Optional[List[SequenceID]] = None,
                 mincomp: Optional[Compatibility] = None,
                 compatibilities: Optional[Dict[SequenceID, Compatibility]] = None,
                 consensus: Optional[SequencePath] = None):
        self.id: AffinityNodeID = id
        self.parent: AffinityNodeID = parent
        self.children: List[AffinityNodeID] = children if children else []
        self.sequences: List[SequenceID] = sequences if sequences else []
        self.mincomp: Compatibility = mincomp if mincomp else Compatibility(0)
        self.compatibilities: Dict[SequenceID, Compatibility] = compatibilities if compatibilities else {}
        self.consensus: SequencePath = consensus

    def __str__(self):
        return f"ID: {self.id}, "\
            f"parent: {self.parent}, " \
            f"children: {self.children}, " \
            f"mincomp: {self.mincomp.value}, " \
            f"consensuss length: {len(self.consensus)}, "\
            f"sequences: {self.sequences}."

    def __eq__(self, other):
        return self.parent == other.parent and \
               self.id == other.id and \
               self.sequences == other.sequences and \
               math.isclose(self.mincomp.valuem, other.mincomp.value) and \
               self.compatibilities == other.compatibilities and \
               self.consensus == other.consensus


class AffinityTree:
    """
    Affinity Tree storage.
    """

    def __init__(self):
        self.nodes: List[AffinityNode] = []

    def get_node(self, node_id: AffinityNodeID) -> AffinityNode:
        """Returns affinity npode with given ID."""

        for node in self.nodes:
            if node.id == node_id:
                return node
        raise AffinityTreeException("No node with given ID.")

    def get_max_node_id(self):

        if len(self.nodes) == 0:
            return -1
        return max([node.id for node in self.nodes])

    def as_newick(self, seq_id_to_seq_name: Dict[SequenceID, str] = None, expand_leaves=False):
        return self._convert_to_newick(seq_id_to_seq_name, expand_leaves)

    def _convert_to_newick(self, seq_id_to_metadata: Dict[SequenceID, str] = None, expand_newick = False) -> str:
        def newick_nhx(newick_tree):

            label = newick_tree.name or ''
            if newick_tree._length:
                for cn in sorted_nodes:
                    if str(cn.id) == newick_tree.name:
                        if seq_id_to_metadata:
                            if len(cn.sequences) == 1:
                                name = seq_id_to_metadata[cn.sequences[0]]["name"]
                                group = seq_id_to_metadata[cn.sequences[0]]["group"]
                                seqid = cn.sequences[0]
                                metadata = f"[&&NHX:name={name}:group={group}:seqid={seqid}:mincomp={cn.mincomp}]"
                            elif len(cn.sequences) == 0:
                                name = f"EmptyAffinityNode {cn.id}"
                                metadata = f"[&&NHX:name={name}:mincomp={cn.mincomp}]"
                            else:
                                name = f"AffinityNode {cn.id}"
                                metadata = f"[&&NHX:name={name}:mincomp={cn.mincomp}]"
                        else:
                            if len(cn.sequences) == 1:
                                name = cn.sequences[0]
                            elif len(cn.sequences) == 0:
                                name = f"EmptyAffinityNode {cn.id}"
                            else:
                                name = f"AffinityNode {cn.id}"
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

        sorted_nodes = sorted(self.nodes, key=lambda x: x.id)

        if expand_newick:
            new_leaves_count = 0
            for node in self.nodes:
                if len(node.children) == 0 and len(node.sequences) > 1:
                    for seq_id in node.sequences:
                        affinity_node_id = len(self.nodes) + new_leaves_count
                        node.children.append(affinity_node_id)
                        sorted_nodes.append(AffinityNode(id=AffinityNodeID(affinity_node_id),
                                                         parent=node.id,
                                                         children=[],
                                                         sequences=[seq_id],
                                                         mincomp=Compatibility(1.0)
                                                         ))
                        new_leaves_count += 1

        nodes_to_process = [(None, sorted_nodes[0])]
        newick_tree = None
        while nodes_to_process:
            n = nodes_to_process.pop()
            node_parent_label = n[0]
            node = n[1]

            label = str(node.id)
            if node.parent is None:
                length = "1"
            else:
                parent_minComp = sorted_nodes[node.parent].mincomp.base_value().value
                length = str((1 - parent_minComp) - (1 - node.mincomp.base_value().value))

            newick_node = Node(name=label, length=length)

            if newick_tree is None:
                newick_tree = newick_node
            else:
                parent_node = newick_tree.get_node(node_parent_label)
                parent_node.add_descendant(newick_node)

            for child in node.children:
                nodes_to_process.append((label, sorted_nodes[child]))
        return "(" + newick_nhx(newick_tree) + ")"
