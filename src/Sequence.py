class Sequence(object):
    def __init__(self, ID, name, title, active = True, nodes_IDs = None):
        self.ID = ID
        self.name = name
        self.title = title
        self.active = active
        self.nodes_IDs = nodes_IDs if nodes_IDs else []

    def __str__(self):
        return """ID: {0},\t name: {1},\t title: {2}, \t active: {3},\t nodes IDs: {4}""".format(
            self.ID,
            self.name,
            self.title,
            self.active,
            self.nodes_IDs)

    def __eq__(self, other):
        return     (self.ID == other.ID
                    and self.name == other.name
                    and self.title == other.title
                    and self.active == other.active
                    and self.nodes_IDs == other.nodes_IDs)

    def add_node_ID(self, node_ID):
        self.nodes_IDs.append(node_ID)

class Source(Sequence):
    def __init__(self, ID, name, title,active = True, nodes_IDs = None, consensusID = -1, weight = -1):
        Sequence.__init__(self, ID=ID, name=name, title=title, active=active, nodes_IDs=nodes_IDs)
        self.consensusID = consensusID
        self.weight = weight

    def __str__(self):
        return Sequence.__str__(self) + """\tconsensusID: {0},\tweight: {1}""".format(
            self.consensusID,
            self.weight)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and self.consensusID == other.consensusID \
                                            and self.weight == other.weight



class Consensus(Sequence):
    def __init__(self, ID, name, title, active = True, nodes_IDs = None, compatibility_to_sources = None):
        Sequence.__init__(self, ID=ID, name=name, title=title, active=active, nodes_IDs=nodes_IDs)
        self.compatibility_to_sources = compatibility_to_sources if compatibility_to_sources else {}

    def __str__(self):
        return Sequence.__str__(self) + """ compatibility_to_sources: {0}""".format(  self.compatibility_to_sources)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and self.compatibility_to_sources == other.compatibility_to_sources