from typing import NewType, Optional

from poapangenome.datamodel.builders.PoagraphBuildException import PoagraphBuildException

NodeID = NewType('NodeID', int)
ColumnID = NewType('ColumnID', int)
BlockID = NewType('BlockID', int)


class Base:
    def __init__(self, b: str):
        if len(b) == 0:
            raise PoagraphBuildException("Empty string for base.")
        if len(b) > 1:
            raise PoagraphBuildException("Base string must have length 1.")
        self.value: bytes = str.encode(b)

    def __eq__(self, other: 'Base'):
        return other and self.value.__eq__(other.value)

    def as_str(self) -> str:
        """Return base value as string."""

        return self.value.decode("ASCII")


class Node:
    def __init__(self,
                 node_id: NodeID,
                 base: Base,
                 aligned_to: Optional[NodeID] = None,
                 column_id: Optional[ColumnID] = None,
                 block_id: Optional[BlockID] = None):
        self.node_id: NodeID = node_id
        self.base: Base = base
        self.aligned_to: NodeID = aligned_to
        self.column_id: ColumnID = column_id
        self.block_id: BlockID = block_id

    def get_base(self):
        return self.base.value.decode("ASCII")

    def __eq__(self, other: 'Node'):
        return (self.node_id == other.node_id
                and self.base == other.base
                and self.aligned_to == other.aligned_to
                # and self.column_id == other.column_id
                # and self.block_id == other.block_id
                )

    def __str__(self):
        return \
            f"id: {self.node_id}, " \
            f"base: {self.get_base()}, " \
            f"aligned_to: {self.aligned_to}, " \
            f"column_id: {self.column_id}, " \
            f"block_id: {self.block_id}"

    def __repr__(self):
        return self.__str__()
