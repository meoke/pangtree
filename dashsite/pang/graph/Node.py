from typing import List, Union, Set
from .nucleotides import decode


class Node:
    def __init__(self,
                 id: int,
                 base: int,
                 in_nodes: Union[Set[int], List[int]],
                 aligned_to: List[int]):
        self.id = id
        self.base = base
        self.in_nodes = in_nodes
        self.aligned_to = aligned_to

    def __eq__(self, other):
        return (self.id == other.id
            and self.base == other.base
            and sorted(self.in_nodes) == sorted(other.in_nodes)
            and self.aligned_to == other.aligned_to)

    def __str__(self):
        return f"id: {self.id}, base: {decode(self.base)}, in_nodes: {self.in_nodes}, aligned_to: {self.aligned_to}"


