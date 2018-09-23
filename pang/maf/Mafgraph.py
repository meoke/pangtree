from .Block import Block
from .Direction import Direction
from . import sorter

class Mafgraph:
    def __init__(self, maf, remove_cycles: bool):
        #todo typ dla mafa
        self.blocks = self._maf_to_blocks(maf)
        if remove_cycles:
            self._sort_blocks()

    def __iter__(self):
        for block in sorted(self.blocks, key=lambda block: block.order_id):
            yield block

    def _maf_to_blocks(self, maf):
        blocks = []
        for i, block in enumerate(maf):
            b = Block(block_id=i, order_id=i, direction=Direction.LEFT_TO_RIGHT, alignment=block)
            blocks.append(b)
        return blocks

    def _sort_blocks(self):
        return sorter.sort_mafblocks(self.blocks)

# pytania:
# czy chce przechowywać treść z AlignIO?

