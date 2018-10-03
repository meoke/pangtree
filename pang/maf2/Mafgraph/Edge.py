from typing import List
from .EdgeType import EdgeType


class Edge(object):
    def __init__(self, to: int, edge_type: EdgeType, sequences=List[int]):
        self.to = to
        self.edge_type = edge_type
        self.sequences = sequences
        self.active = True

    def __str__(self):
        return f"Edge to: {self.to}, type: {self.edge_type}, sequences: {[str(s.seq) for s in self.sequences]}"
