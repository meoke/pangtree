from typing import List
from . import PathManager
from .Node import Node


NodesList = List[Node]


class Graph:
    def __init__(self, nodes_count: int=0, start_node_id: int = 0, nodes: NodesList=None):
        if nodes:
            self.nodes = nodes
        else:
            # todo metoda do sprawdzania poprawności
            self.nodes = [None]*nodes_count

    def __eq__(self, other):
        return self.nodes == other.nodes

    def update_nodes(self, new_nodes: NodesList):
        if new_nodes:
            self.nodes[new_nodes[0].id: new_nodes[-1].id] = new_nodes

    def trim(self, nodes_count: int):
        del self.nodes[nodes_count:]

