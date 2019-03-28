from pangraph.custom_types import NodeID, Base, BlockID, ColumnID


class Node:
    def __init__(self,
                 node_id: NodeID,
                 base: Base,
                 aligned_to: NodeID,
                 column_id: ColumnID = None,
                 block_id: BlockID = None):
        self.id = node_id
        self.base = base
        self.aligned_to = aligned_to
        self.column_id = column_id
        self.block_id = block_id

    def get_base(self):
        return self.base.decode("ASCII")

    def __eq__(self, other):
        return (self.id == other.id
                and self.base == other.base
                and self.aligned_to == other.aligned_to)

    def __str__(self):
        base = self.base.decode("ASCII")
        return \
            f"id: {self.id}, " \
            f"base: {base}, " \
            f"aligned_to: {self.aligned_to}, " \
            f"column_id: {self.column_id}, " \
            f"block_id: {self.block_id}"



