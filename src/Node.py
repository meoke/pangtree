class Node(object):
    def __init__(self, ID, base, in_nodes=None, aligned_to=None, sources_count = 0, consensuses_count = 0):
        self.ID = ID
        self.base = base
        self.in_nodes = in_nodes if in_nodes else set()
        self.aligned_to = aligned_to if aligned_to else set()
        self.sources_count = sources_count
        self.consensuses_count = consensuses_count

    def __eq__(self, other):
        return (self.ID == other.ID
            and self.base == other.base
            and self.in_nodes == other.in_nodes
            and self.aligned_to == other.aligned_to
            and self.sources_count == other.sources_count
            and self.consensuses_count == other.consensuses_count)


    def __str__(self):
        return  """ID: {0},\t base: {1},\t in_nodes: {2},\t aligned_to: {3},\t sources_count: {4},\t consensuses_count: {5}""".format(
            self.ID,
            self.base,
            self.in_nodes,
            self.aligned_to,
            self.sources_count,
            self.consensuses_count)
