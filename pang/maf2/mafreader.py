from .Mafgraph.Block import Block
from .Mafgraph.EdgeType import EdgeType


def read_mafblocks(maf):
    blocks = []
    for i, mafblock in enumerate(maf):
        b = Block(block_id=i, order_id=i, alignment=mafblock, reversed=False)
        for sequence in mafblock:
            next_block_id, edge_type = _find_next_block_id(maf, sequence)
            b.add_out_edge(to=next_block_id, sequence=sequence, edge_type=edge_type)
        blocks.append(b)
    return blocks


def _find_next_block_id(maf, sequence):
    return -1, EdgeType.RIGHT_TO_LEFT
