from typing import List

from datamodel.input_types import Maf
from mafgraph.sorter import sort_mafblocks


class DAGMafNode:
    def __init__(self, block_id, alignment, orient, order, out_edges):
        self.id = block_id
        self.alignment = alignment
        self.orient = orient  # orientation relative to parent
        self.order = order
        self.out_edges = out_edges


class DAGMaf:
    def __init__(self, dagmaf_nodes: List[DAGMafNode]):
        self.dagmaf_nodes = dagmaf_nodes


def get_dagmaf(maf: Maf) -> DAGMaf:
    sorted_blocks = sort_mafblocks(maf.filecontent)
    dagmafnodes = [
        DAGMafNode(block_id=b.id,
                   alignment=b.alignment,
                   orient=b.orient,
                   order=b.order(),
                   out_edges=b.out_edges)
        for b in sorted_blocks
    ]
    return DAGMaf(dagmafnodes)
