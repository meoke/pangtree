from typing import List
from .EdgeType import EdgeType
from .SequenceInfo import SequenceInfo


class Edge(object):
    def __init__(self, to: int, edge_type: EdgeType, sequences=List[SequenceInfo]):
        self.to = to
        self.edge_type = edge_type
        self.sequences = sequences

    def __str__(self):
        return f"Edge to: {self.to}, type: {self.edge_type}, sequences: {[str(s) for s in self.sequences]}"
