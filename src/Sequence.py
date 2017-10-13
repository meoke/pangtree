class Sequence(object):
    def __init__(self, ID, name, title, first_node_ID = -1, active = True, nodes_IDs = None):
        self.ID = ID
        self.name = name
        self.title = title
        self.first_node_id = first_node_ID
        self.active = active
        self.nodes_IDs = nodes_IDs if nodes_IDs else set()

    def __str__(self):
        return """ID: {0},\t name: {1},\t title: {2},\t first_node_id: {3},\t active: {4},\t nodes IDs: {5}""".format(
            self.ID,
            self.name,
            self.title,
            self.first_node_id,
            self.active,
            self.nodes_IDs)

    def __eq__(self, other):
        return     (self.ID == other.ID
                    and self.name == other.name
                    and self.title == other.title
                    and self.first_node_id == other.first_node_id
                    and self.active == other.active
                    and self.nodes_IDs == other.nodes_IDs)

    def add_node_ID(self, node_ID):
        self.nodes_IDs.add(node_ID)

class Source(Sequence):
    def __init__(self, ID, name, title, first_node_ID = -1, active = True, nodes_IDs = set(), consensusID = -1, weight = -1):
        Sequence.__init__(self, ID=ID, name=name, title=title, first_node_ID=first_node_ID, active=active, nodes_IDs=nodes_IDs)
        self.consensusID = consensusID
        self.weight = weight

    def __str__(self):
        return Sequence.__str__(self) + """\tconsensusID: {0},\tweight: {1}""".format(
            self.consensusID,
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