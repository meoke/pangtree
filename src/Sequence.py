class Sequence(object):
    def __init__(self, ID, name, title, nodes_count = 0, first_node_ID = -1, active = True, nodes_IDs = set()):
        self.ID = ID
        self.name = name
        self.title = title
        self.nodes_count = nodes_count
        self.first_node_id = first_node_ID
        self.active = active
        self.nodes_IDs = nodes_IDs

    def __str__(self):
        return """ID: {0},
                name: {1},
                title: {2},
                nodes_count: {3},
                first_node_id: {4},
                active: {5},
                nodes IDs: {6}""".format(self.ID,
                                         self.name,
                                         self.title,
                                         self.nodes_count,
                                         self.first_node_id,
                                         self.active,
                                         self.nodes_IDs)

    def __eq__(self, other):
        return     (self.ID == other.ID
                    and self.name == other.name
                    and self.title == other.title
                    and self.nodes_count == other.nodes_count
                    and self.first_node_id == other.first_node_id
                    and self.active == other.active
                    and self.nodes_IDs == other.nodes_IDs)

    def add_node_ID(self, node_ID):
        self.nodes_IDs.add(node_ID)

class Source(Sequence):
    def __init__(self, ID, name, title, nodes_count = 0, first_node_ID = -1, active = True, nodes_IDs = set(), consensusID = -1, weight = -1):
        Sequence.__init__(self, ID, name, nodes_count, first_node_ID, title, active, nodes_IDs)
        self.consensusID = consensusID
        self.weight = weight

    def __str__(self):
        return Sequence.__str__(self) + """ consensusID: {0},
                                            weight: {1}""".format(  self.consensusID,
                                                                    self.weight)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and self.consensusID == other.consensusID \
                                            and self.weight == other.weight


# class Consensus(Sequence):
#     def __init__(self, ID, name, nodes_count = 0, first_node_ID = -1, title = "", active = True, compatibility_to_sources = {}):
#         Sequence.__init__(self, ID, name, nodes_count, first_node_ID, title, active)
#         self.compatibility_to_sources = compatibility_to_sources
#
#     def __str__(self):
#         return Sequence.__str__(self) + """ compatibility_to_sources: {0}""".format(  self.compatibility_to_sources)
#
#     def __eq__(self, other):
#         return Sequence.__eq__(self, other) and self.compatibility_to_sources == other.compatibility_to_sources