class Node(object):
    def __init__(self, ID, base, in_nodes=set(), aligned_to=set(), sources_count = 0, consensuses_count = 0):
        self.ID = ID
        self.base = base
        self.in_nodes = in_nodes
        self.aligned_to = aligned_to
        self.sources_count = sources_count
        self.consensuses_count = consensuses_count

    def __eq__(self, other):
        return (self.ID == other.ID
            and self.base == other.base
            and self.in_nodes == other.in_nodes
            and self.aligned_to == other.alignedTo)


    def __str__(self):
        return  """ID: {0},
                   base: {1},
                   in_nodes: {2},
                   aligned_to: {3}""".format(self.ID,
                                             self.base,
                                             self.in_nodes,
                                             self.aligned_to)
