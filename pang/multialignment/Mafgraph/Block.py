from .Edge import Edge
from .EdgeType import EdgeType


class Block:
    def __init__(self, block_id: int, order_id: int, alignment):
        self.block_id = block_id
        self.order_id = order_id
        self.alignment = alignment
        self.out_edges = []

    def __str__(self):
        block_str = [f"Block id: {self.block_id}, order id: {self.order_id}",
                     f"Out edges: {[str(e) for e in self.out_edges]}"]
        return "\n".join(block_str)

    def add_out_edge(self, to: int, sequence, edge_type: EdgeType):
        for edge in self.out_edges:
            if edge.to == to and edge.edge_type == edge_type:
                edge.sequences.append(sequence)
                return
        self.out_edges.append(Edge(to=to, edge_type=edge_type, sequences=[sequence]))
