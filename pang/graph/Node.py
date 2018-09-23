from typing import List


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

