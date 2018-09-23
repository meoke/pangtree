from typing import List
from . import PathManager
from .Node import Node


NodesList = List[Node]


class Graph:
    def __init__(self, nodes_count):
        self.nodes = [None]*nodes_count

    def update_nodes(self, new_nodes: NodesList):
        if new_nodes:
            self.nodes[new_nodes[0].id: new_nodes[-1].id] = new_nodes

