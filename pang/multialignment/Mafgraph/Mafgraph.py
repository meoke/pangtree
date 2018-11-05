from maf2 import sorter
from maf2 import mafreader


class Mafgraph(object):
    def __init__(self, maf, remove_cycles: bool):
        self.blocks = self._maf_to_blocks(maf)
        if remove_cycles:
            self._sort_blocks()

    def __iter__(self):
        for block in sorted(self.blocks, key=lambda b: b.order_id):
            yield block

    def __str__(self):
        mafgraph_str = []
        for i, block in enumerate(self.blocks):
            mafgraph_str.append(f"Block {i}:")
            mafgraph_str.append(f"{block}")
        return "\n".join(mafgraph_str)

    def _maf_to_blocks(self, maf):
        return mafreader.read_mafblocks(maf)

    def _sort_blocks(self):
        return sorter.sort_mafblocks(self.blocks)


