from mafgraph.sorter import sort_mafblocks
from mafgraph.mafreader import read_maf

from mafgraph.graph.Block import Block


class DAGMafNode:
    def __init__(self, block_id, alignment, orient, order, out_edges):
        self.id = block_id
        self.alignment = alignment
        self.orient = orient  # orientation relative to parent
        self.order = order
        self.out_edges = out_edges


class DAGMaf:
    def __init__(self, multialignment):
        sorted_blocks = sort_mafblocks(multialignment)
        self.dagmafnodes = [
            DAGMafNode(block_id=b.id,
                       alignment=b.alignment,
                       orient=b.orient,
                       order=b.order(),
                       out_edges=b.out_edges)
            for b in sorted_blocks
        ]
