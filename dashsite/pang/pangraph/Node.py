from typing import List, Union, Set
from .nucleotides import decode


class Node:
    def __init__(self,
                 node_id: int,
                 base: int,
                 aligned_to: List[int],
                 column_id: int = None,
                 block_id: int = None):
        self.id = node_id
        self.base = base
        self.aligned_to = aligned_to
        self.column_id = column_id
        self.block_id = block_id

    def __eq__(self, other):
        return (self.id == other.id
                and self.base == other.base
                and self.aligned_to == other.aligned_to)

    def __str__(self):
        return \
            f"id: {self.id}, " \
            f"base: {decode(self.base)}, " \
            f"aligned_to: {self.aligned_to}, " \
            f"column_id: {self.column_id}, " \
            f"block_id: {self.block_id}"


