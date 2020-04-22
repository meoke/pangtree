"""AffinityTree and related objects definitions."""

import math
from typing import List, NewType, Dict, Optional

import newick
from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import msa

AffinityNodeID = NewType('AffinityNodeID', int)


class AffinityNode(object):
    """Node of Affinity Tree.

    Args:
        id_: Affinity Node ID.
        parent: ID of the parent node.
        children: IDs of the children nodes.
        sequences: IDs of the _sequences assigned to this node.
        mincomp: Minimum from the compatibilities of this consensus to the
            assigned sequences
        compatibilities: Dictionary of compatibilities to all sequences.
        consensus: Path of the consensus defined as path in Poagraph.

    Attributes:
        id_ (AffinityNodeID): ID of this node
        parent (AffinityNodeID): ID of the parent node
        children(List[AffinityNodeID]): IDs of the children nodes.
        sequences(List[SequenceID]): IDs of the sequences assigned to the node.
        mincomp(Compatibility): Minimum from the compatibilities of this
            consensus to the assigned sequences
        compatibilities(Dict[SequenceID, Compatibility]): Dictionary of
            compatibilities to any _sequences.
        consensus(SeqPath): Path of the consensus defined as path in Poagraph.
    """

    def __init__(self,
                 id_: AffinityNodeID,
                 parent: Optional[AffinityNodeID] = None,
                 children: Optional[List[AffinityNodeID]] = None,
                 sequences: Optional[List[msa.SequenceID]] = None,
                 mincomp: Optional[graph.Compatibility] = None,
                 compatibilities: Optional[Dict[msa.SequenceID,
                                                graph.Compatibility]] = None,
                 consensus: Optional[graph.SeqPath] = None):
        self.id_: AffinityNodeID = id_
        self.parent: AffinityNodeID = parent
        self.children: List[AffinityNodeID] = children if children else []
        self.sequences: List[msa.SequenceID] = sequences if sequences else []
        self.mincomp: graph.Compatibility = mincomp if mincomp else graph.Compatibility(0)
        self.compatibilities: Dict[msa.SequenceID,
                                   graph.Compatibility] = compatibilities if compatibilities else {}
        self.consensus: graph.SeqPath = consensus

    def __str__(self):
        return f"ID: {self.id_}, "\
               f"parent: {self.parent}, " \
               f"children: {self.children}, " \
               f"mincomp: {str(self.mincomp)}, " \
               f"consensus length: {len(self.consensus)}, "\
               f"_sequences: {self.sequences}."

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
        nodes (List[AffinityNode]): All nodes of the tree. Relations between
            nodes are described by their attributes.
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
            The largest AffinityNodeID used in this Affinity Tree. Returns
                AffinityNodeID(-1) if the tree has no nodes.
        """

        if len(self.nodes) == 0:
            return AffinityNodeID(-1)
        return max([node.id_ for node in self.nodes])

    def as_newick(self,
                  seq_id_to_metadata: Dict[msa.SequenceID,
                                           graph.SequenceMetadata] = None,
                  separate_leaves=False) -> str:
        """Returns Affinity Tree in Newick format.

        Args:
            seq_id_to_metadata: Dictionary of _sequences IDs to the desired
                name used in newick file. For example:
                                {SequenceID('KM0123'): 'cat',
                                SequenceID('ZX124'): 'dog'}
            separate_leaves: A switch to control if tree leaves having
                assigned multiple _sequences should have appended
                children singleton leaves single sequence assigned.

        Returns:
            A string with the Affinity Tree converted to newick format.
            https://en.wikipedia.org/wiki/Newick_format
            If the tree has no nodes, an empty string is returned.
        """

        def _get_sequence_attr_if_exists(seq_metadata: graph.SequenceMetadata,
                                         attr: str) -> str:
            """Returns dictionary value if they key attr exists."""

            if attr in seq_metadata:
                return str(seq_metadata[attr])
            else:
                return ""

        def _newick_nhx(newick_node: newick.Node) -> str:
            """Converts newick tree to newick string"""

            node_label = newick_node.name or ''
            if newick_node._length:
                for cn in sorted_nodes:
                    if str(cn.id_) == newick_node.name:
                        if seq_id_to_metadata:
                            if len(cn.sequences) == 1:
                                name = _get_sequence_attr_if_exists(seq_id_to_metadata[cn.sequences[0]], "name")
                                if name == "":
                                    name = cn.sequences[0]
                                group = _get_sequence_attr_if_exists(seq_id_to_metadata[cn.sequences[0]], "group")
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
                except Exception:
                    print("metadata")
            descendants = ','.join([_newick_nhx(n)
                                    for n in newick_node.descendants])
            if descendants:
                descendants = '(' + descendants + ')'
            return descendants + node_label

        if not self.nodes:
            return ""

        sorted_nodes = sorted(self.nodes, key=lambda x: x.id_)
        remove_children = []
        if separate_leaves:
            new_leaves_count = 0
            for node in self.nodes:
                if len(node.children) == 0 and len(node.sequences) > 1:
                    for seq_id in node.sequences:
                        affinity_node_id = len(self.nodes) + new_leaves_count
                        node.children.append(affinity_node_id)
                        remove_children.append(node.id_)
                        sorted_nodes.append(AffinityNode(id_=AffinityNodeID(affinity_node_id),
                                                         parent=node.id_,
                                                         children=[],
                                                         sequences=[seq_id],
                                                         mincomp=graph.Compatibility(1.0)
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

            newick_node = newick.Node(name=label, length=length)

            if newick_tree is None:
                newick_tree = newick_node
            else:
                parent_node = newick_tree.get_node(node_parent_label)
                parent_node.add_descendant(newick_node)

            for child in node.children:
                nodes_to_process.append((label, sorted_nodes[child]))
        for node in self.nodes:
            if node.id_ in remove_children:
                node.children = []

        return "(" + _newick_nhx(newick_tree) + ")"
