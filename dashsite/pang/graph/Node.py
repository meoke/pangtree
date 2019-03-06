from typing import List, Union, Set
from .nucleotides import decode


class Node:
    def __init__(self,
                 id: int,
                 base: int,
                 in_nodes: Union[Set[int], List[int]],
                 aligned_to: List[int],
                 column_id: int = None,
                 block_id: int = None):
        self.id = id
        self.base = base
        self.in_nodes = in_nodes
        self.aligned_to = aligned_to
        self.column_id = column_id
        self.block_id = block_id

    def __eq__(self, other):
        return (self.id == other.id
            and self.base == other.base
            and sorted(self.in_nodes) == sorted(other.in_nodes)
            and self.aligned_to == other.aligned_to)

    def __str__(self):
        return \
            f"id: {self.id}, " \
            f"base: {decode(self.base)}, " \
            f"in_nodes: {self.in_nodes}, " \
            f"aligned_to: {self.aligned_to}, " \
            f"column_id: {self.column_id}, " \
            f"block_id: {self.block_id}"


