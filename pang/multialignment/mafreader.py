from .Mafgraph.Block import Block
from .Mafgraph.EdgeType import EdgeType
from .Mafgraph.SequenceInfo import SequenceInfo


def read_mafblocks(maf):
    """Read blocks from maf file."""
    blocks = []
    for i, mafblock in enumerate(maf):
        b = Block(block_id=i, order_id=i, alignment=mafblock)
        for sequence in mafblock:
            next_block_id, edge_type = _find_next_block_id(maf, sequence)
            b.add_out_edge(to=next_block_id, sequence=SequenceInfo(sequence.id, sequence.annotations['start']), edge_type=edge_type)
        blocks.append(b)
    return blocks


def _find_next_block_id(maf, sequence):
    """Find block connected with given sequence"""
    return -1, EdgeType.RIGHT_TO_LEFT
