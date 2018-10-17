from typing import List
from .nucleotides import decode


NodesIDsList = List[int]


class Node:
    def __init__(self,
                 id: int,
                 base: int,
                 in_nodes: NodesIDsList,
                 aligned_to: NodesIDsList):
        self.id = id
        self.base = base
        self.in_nodes = in_nodes
        self.aligned_to = aligned_to

    def __eq__(self, other):
        return (self.id == other.id
            and self.base == other.base
            and self.in_nodes == other.in_nodes
            and self.aligned_to == other.aligned_to)

    def __str__(self):
        return f"id: {self.id}, base: {decode(self.base)}, in_nodes: {self.in_nodes}, aligned_to: {self.aligned_to}"


