import numpy as np
from bisect import insort_left

class Node(object):
    def __init__(self, ID, base, in_nodes=np.array([]), aligned_to=None, consensuses_count=0):
        self.ID = ID
        self.base = base
        self.in_nodes = in_nodes if in_nodes.size else np.array([], dtype=int)
        self.aligned_to = aligned_to
        self.consensuses_count = consensuses_count #TODO czy to potrzebne? sprawdzic wydajnosc bez tego

    def __eq__(self, other):
        return (self.ID == other.ID
                and self.base == other.base
                and np.array_equal(self.in_nodes, other.in_nodes)
                and self.aligned_to == other.aligned_to
                and self.consensuses_count == other.consensuses_count)


    def __str__(self):
        return  """ID: {0},\t base: {1},\t in_nodes: {2},\t aligned_to: {3},\t consensuses_count: {4}""".format(
            self.ID,
            self.base,
            self.in_nodes,
            self.aligned_to,
            self.consensuses_count)

    def add_in_node(self, node_ID):
        # insort_left(self.in_nodes, node_ID)

        # self.in_nodes[idx] = node_ID
        self.in_nodes = np.unique(np.append(self.in_nodes, node_ID))
