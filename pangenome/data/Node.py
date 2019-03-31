from typing import NewType

from data.builders import PoagraphBuildException

NodeID = NewType('NodeID', int)
ColumnID = NewType('ColumnID', int)
BlockID = NewType('BlockID', int)


class Base:
    def __init__(self, b: str):
        if len(b) == 0:
            raise PoagraphBuildException("Empty string for base.")
        if len(b) > 1:
            raise PoagraphBuildException("Base string must have length 1.")
        self.value = str.encode(b)


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
        return self.base.value.decode("ASCII")

    def __eq__(self, other: 'Node'):
        return (self.id == other.id
                and self.base == other.base
                and self.aligned_to == other.aligned_to
                and self.column_id == other.column_id
                and self.block_id == other.block_id)

    def __str__(self):
        return \
            f"id: {self.id}, " \
            f"base: {self.get_base()}, " \
            f"aligned_to: {self.aligned_to}, " \
            f"column_id: {self.column_id}, " \
            f"block_id: {self.block_id}"
