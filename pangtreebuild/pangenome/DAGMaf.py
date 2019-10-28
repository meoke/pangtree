from typing import List


class DAGMafNode:
    """Single node of Multialignment DAG."""

    def __init__(self, block_id, alignment, orient, order, out_edges):
        self.id = block_id
        self.alignment = alignment
        self.orient = orient  # orientation relative to parent
        self.order = order
        self.out_edges = out_edges


class DAGMaf:
    """MAF multialignment converted to DAG."""

    def __init__(self, dagmaf_nodes: List[DAGMafNode]):
        self.dagmaf_nodes = dagmaf_nodes
