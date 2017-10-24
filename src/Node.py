class Node(object):
    def __init__(self, currentID, base, in_nodes=None, aligned_to=None, consensuses_count=0, sources=None):
        self.currentID = currentID
        self.base = base
        self.in_nodes = in_nodes if in_nodes else set() # always global node ID
        self.aligned_to = aligned_to if aligned_to else set() # alwys global node ID
        self.consensuses_count = consensuses_count
        self.sources = sources if sources else set() # always global source ID!

    def __eq__(self, other):
        return (self.currentID == other.currentID
            and self.base == other.base
            and self.in_nodes == other.in_nodes
            and self.aligned_to == other.aligned_to
            and self.sources == other.sources
            and self.consensuses_count == other.consensuses_count)


    def __str__(self):
        return  """currentID: {0},\t base: {1},\t in_nodes: {2},\t aligned_to: {3},\t sources: {4},\t consensuses_count: {5}""".format(
            self.currentID,
            self.base,
            self.in_nodes,
            self.aligned_to,
            self.sources,
            self.consensuses_count)
