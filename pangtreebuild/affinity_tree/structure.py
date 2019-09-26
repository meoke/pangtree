"""AffinityTree and related objects definitions: AffinityNode, AffinityNodeID, Compatibility."""

import math
from typing import List, NewType, Dict, Optional, Union, Any

from pangtreebuild.affinity_tree import parameters
from newick import Node
from pangtreebuild.datamodel import Sequence


AffinityNodeID = NewType('AffinityNodeID', int)


class Compatibility(object):
    """Asymetric similiarity measure of two poagraph paths.

    Args:
        compatibility: Raw compatibility value - count of common nodes devided by length of one of the paths.
        p: Parameter to control compatibility value interpretation. Compatibility is raised to the power of P.

    Attributes:
        value (float): Compatibility value raised to the power of p.
        p (float): P parameter value.

    """

    def __init__(self, compatibility: float, p: parameters.P = parameters.P(1)):
        self.value: float = compatibility**p.value
        self.p: float = p.value

    def _check_p_equality(self, other: Union["Compatibility", Any]) -> None:
        if isinstance(other, Compatibility):
            assert self.p == other.p, 'Cannot compare compatibilities with different p values.'
        else:
            return

    def __eq__(self, other: Union["Compatibility", parameters.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value == other.value

    def __lt__(self, other: Union["Compatibility", parameters.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value < other.value

    def __le__(self, other: Union["Compatibility", parameters.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value <= other.value

    def __gt__(self, other: Union["Compatibility", parameters.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value > other.value

    def __ge__(self, other: Union["Compatibility", parameters.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value >= other.value

    def __sub__(self, other: Union["Compatibility", parameters.Stop]) -> "Compatibility":
        self._check_p_equality(other)
        return Compatibility(self.value - other.value, parameters.P(self.p))

    def __str__(self) -> str:
        return f"""{self.value}"""

    def __repr__(self) -> str:
        return f"""value: {self.value}, p: {self.p}"""

    def base_value(self) -> "Compatibility":
        """Get compatibility value without P transformation.

        Returns:
        Compatibility object with the original compatibility value and P=1."""

        return Compatibility(self.value ** (1 / self.p))


class AffinityNode(object):
    """Node of Affinity Tree.

    Args:
        id_: Affinity Node ID.
        parent: ID of the parent node.
        children: IDs of the children nodes.
        sequences: IDs of the sequences assigned to this node.
        mincomp: Minimum from the compatibilities of this consensus to the assigned sequences
        compatibilities: Dictionary of compatibilities SequenceID:Compatibility to any sequences.
        consensus: Path of the consensus defined as path in corresponding Poagraph.

    Attributes:
        id_ (AffinityNodeID): ID of this node
        parent (AffinityNodeID): ID of the parent node
        children(List[AffinityNodeID]): IDs of the children nodes.
        sequences(List[SequenceID]): IDs of the sequences assigned to this node.
        mincomp(Compatibility): Minimum from the compatibilities of this consensus to the assigned sequences
        compatibilities(Dict[SequenceID, Compatibility]): Dictionary of compatibilities to any sequences.
        consensus(SeqPath): Path of the consensus defined as path in corresponding Poagraph.
    """

    def __init__(self,
                 id_: AffinityNodeID,
                 parent: Optional[AffinityNodeID] = None,
                 children: Optional[List[AffinityNodeID]] = None,
                 sequences: Optional[List[Sequence.SequenceID]] = None,
                 mincomp: Optional[Compatibility] = None,
                 compatibilities: Optional[Dict[Sequence.SequenceID, Compatibility]] = None,
                 consensus: Optional[Sequence.SeqPath] = None):
        self.id_: AffinityNodeID = id_
        self.parent: AffinityNodeID = parent
        self.children: List[AffinityNodeID] = children if children else []
        self.sequences: List[Sequence.SequenceID] = sequences if sequences else []
        self.mincomp: Compatibility = mincomp if mincomp else Compatibility(0)
        self.compatibilities: Dict[Sequence.SequenceID, Compatibility] = compatibilities if compatibilities else {}
        self.consensus: Sequence.SeqPath = consensus

    def __str__(self):
        return f"ID: {self.id_}, "\
               f"parent: {self.parent}, " \
               f"children: {self.children}, " \
               f"mincomp: {str(self.mincomp)}, " \
               f"consensus length: {len(self.consensus)}, "\
               f"sequences: {self.sequences}."

    def __eq__(self, other):
        return self.parent == other.parent and \
               self.id_ == other.id and \
               self.sequences == other.sequences and \
               math.isclose(self.mincomp.value, other.mincomp.value) and \
               self.compatibilities == other.compatibilities and \
               self.consensus == other.consensus


class AffinityTree(object):
    """
    Affinity Tree defined as list of nodes.

    No correctness checking is performed.

    Args:
        nodes: Nodes of the tree.

    Attributes:
        nodes (List[AffinityNode]): All nodes of the tree. Relations between nodes are described by their attributes.
    """

    def __init__(self, nodes: Optional[List[AffinityNode]] = None):
        self.nodes: List[AffinityNode] = nodes if nodes else []

    def get_node(self, id_: AffinityNodeID) -> AffinityNode:
        """Returns affinity node with given ID.

        Args:
            id_: ID of the node to return.

        Returns:
            Affinity Node of given ID.

        Raises:
            KeyError: No node with given ID exists in this Affinity Tree.
        """

        for node in self.nodes:
            if node.id_ == id_:
                return node
        raise KeyError("No node with given ID.")

    def get_max_node_id(self) -> AffinityNodeID:
        """Returns the largest number used as the tree node ID.

        Returns:
            The largest AffinityNodeID used in this Affinity Tree. Returns AffinityNodeID(-1) if the tree has no nodes.
        """

        if len(self.nodes) == 0:
            return AffinityNodeID(-1)
        return max([node.id_ for node in self.nodes])

    def as_newick(self, seq_id_to_metadata: Dict[Sequence.SequenceID, Sequence.SequenceMetadata] = None, separate_leaves=False) -> str:
        """Returns Affinity Tree in Newick format.

        Args:
            seq_id_to_metadata: Dictionary of sequences IDs to the desired name used in newick file. For example:
                                {SequenceID('KM0123'): 'cat',
                                SequenceID('ZX124'): 'dog'}
            separate_leaves: A switch to control if tree leaves having assigned multiple sequences should have appended
                             children singleton leaves single sequence assigned.

        Returns:
            A string with the Affinity Tree converted to newick format. https://en.wikipedia.org/wiki/Newick_format
            If the tree has no nodes, an empty string is returned.
        """

        def _newick_nhx(newick_node: Node) -> str:
            """Converts newick tree to newick string"""

            node_label = newick_node.name or ''
            if newick_node._length:
                for cn in sorted_nodes:
                    if str(cn.id_) == newick_node.name:
                        if seq_id_to_metadata:
                            if len(cn.sequences) == 1:
                                name = seq_id_to_metadata[cn.sequences[0]]["name"]
                                group = seq_id_to_metadata[cn.sequences[0]]["group"]
                                seqid = cn.sequences[0]
                                metadata = f"[&&NHX:name={name}:group={group}:seqid={seqid}:mincomp={cn.mincomp}]"
                            elif len(cn.sequences) == 0:
                                name = f"EmptyAffinityNode {cn.id_}"
                                metadata = f"[&&NHX:name={name}:mincomp={cn.mincomp}]"
                            else:
                                name = f"AffinityNode {cn.id_}"
                                metadata = f"[&&NHX:name={name}:mincomp={cn.mincomp}]"
                        else:
                            if len(cn.sequences) == 1:
                                name = cn.sequences[0]
                            elif len(cn.sequences) == 0:
                                name = f"EmptyAffinityNode {cn.id_}"
                            else:
                                name = f"AffinityNode {cn.id_}"
                            mincomp = cn.mincomp
                            metadata = f"[&&NHX:name={name}:mincomp={mincomp}]"
                try:
                    node_label += ':' + newick_node._length + metadata
                except:
                    print("metadata")
            descendants = ','.join([_newick_nhx(n) for n in newick_node.descendants])
            if descendants:
                descendants = '(' + descendants + ')'
            return descendants + node_label

        if not self.nodes:
            return ""

        sorted_nodes = sorted(self.nodes, key=lambda x: x.id_)

        if separate_leaves:
            new_leaves_count = 0
            for node in self.nodes:
                if len(node.children) == 0 and len(node.sequences) > 1:
                    for seq_id in node.sequences:
                        affinity_node_id = len(self.nodes) + new_leaves_count
                        node.children.append(affinity_node_id)
                        sorted_nodes.append(AffinityNode(id_=AffinityNodeID(affinity_node_id),
                                                         parent=node.id_,
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

            label = str(node.id_)
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
        return "(" + _newick_nhx(newick_tree) + ")"
