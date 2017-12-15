class Sequence(object):
    def __init__(self, currentID, name, title, active = True, nodes_IDs = None):
        self.currentID = currentID
        self.name = name
        self.title = title
        self.active = active
        self.nodes_IDs = nodes_IDs if nodes_IDs else []

    def __str__(self):
        return """ID: {0},\t name: {1},\t title: {2}, \t active: {3},\t nodes IDs: {4}""".format(
            self.currentID,
            self.name,
            self.title,
            self.active,
            self.nodes_IDs)

    def __eq__(self, other):
        return     (self.currentID == other.currentID
                    and self.name == other.name
                    and self.title == other.title
                    and self.active == other.active
                    and self.nodes_IDs == other.nodes_IDs)

    def add_node_ID(self, node_globalID):
        self.nodes_IDs.append(node_globalID)

class Source(Sequence):
    def __init__(self, currentID, name, title, active = True, nodes_IDs = None, consensusID = -1, weight = -1):
        Sequence.__init__(self, currentID=currentID, name=name, title=title, active=active, nodes_IDs=nodes_IDs)
        self.consensusID = consensusID
        self.consensuses = []
        self.weight = weight

    def __str__(self):
        return Sequence.__str__(self) + """\tconsensusID: {0},\tconsensuses: {1},\tweight: {2}""".format(
            self.consensusID,
            self.consensuses,
            self.weight)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and self.consensusID == other.consensusID \
                                            and self.consensuses == other.consensuses \
                                            and self.weight == other.weight



class Consensus(Sequence):
    def __init__(self, currentID, name, title, active=True, nodes_IDs=None, compatibility_to_sources=None, sources_IDs=None):
        Sequence.__init__(self, currentID=currentID, name=name, title=title, active=active, nodes_IDs=nodes_IDs)
        self.compatibility_to_sources = compatibility_to_sources if compatibility_to_sources else {}
        self.level = -1
        self.sources_IDs = sources_IDs if sources_IDs else []
        self.parent_consensus = []

    def __str__(self):
        return Sequence.__str__(self) + """ compatibility_to_sources: {0}""".format(  self.compatibility_to_sources)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and self.compatibility_to_sources == other.compatibility_to_sources