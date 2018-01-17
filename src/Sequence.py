import numpy as np

class Sequence(object):
    def __init__(self, ID, name, title, nodes_IDs = np.array([]), max_nodes_count=0):
        self.ID = ID
        self.name = name
        self.title = title
        if nodes_IDs.size:
            self.nodes_IDs = nodes_IDs
            self.nodes_count = nodes_IDs.size
        else:
            self.nodes_IDs = np.zeros((max_nodes_count),dtype=np.uint32) #nodes_IDs #TODO czy to na pewno sie nie krzaczy
            self.nodes_count = 0

    def __str__(self):
        return """ID: {0},\t name: {1},\t title: {2}, \t nodes IDs: {3}""".format(
            self.ID,
            self.name,
            self.title,
            self.nodes_IDs)

    def __eq__(self, other):
        return     (self.ID == other.ID
                    and self.name == other.name
                    and self.title == other.title
                    and np.array_equal(self.nodes_IDs, other.nodes_IDs))

    def truncate_nodes_ids(self):
        self.nodes_IDs.resize(self.nodes_count)
    # def add_node(self, node_ID):
    #     self.nodes_IDs = np.append(self.nodes_IDs, node_ID)


class Source(Sequence):
    def __init__(self, ID, name, title, nodes_IDs = np.array([]), weight = -1, max_nodes_count=0):
        Sequence.__init__(self, ID=ID, name=name, title=title, nodes_IDs=nodes_IDs, max_nodes_count=max_nodes_count)
        #self.consensuses = {}
        self.weight = weight

    def __str__(self):
        return Sequence.__str__(self) + """\tweight: {0}""".format(
            self.weight)
        # todo jeśli nie używamy consensuses w source, to wystarczy powyższe
        # return Sequence.__str__(self) + """\tconsensusID: {0},\tconsensuses: {1},\tweight: {2}""".format(
        #     self.consensusID,
        #     self.consensuses,
        #     self.weight)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and self.weight == other.weight
        # todo jeśli nie używamy consensuses w source, to wystarczy powyższe
        # return Sequence.__eq__(self, other) and self.consensusID == other.consensusID \
        #                                     and self.consensuses == other.consensuses \
        #                                     and self.weight == other.weight



# class Consensus(Sequence):
#     def __init__(self, currentID, name, title, active=True, nodes_IDs=None, compatibility_to_sources=None, sources_IDs=None):
#         Sequence.__init__(self, currentID=currentID, name=name, title=title, active=active, nodes_IDs=nodes_IDs)
#         self.compatibility_to_sources = compatibility_to_sources if compatibility_to_sources else {}
#         self.level = -1
#         self.sources_IDs = sources_IDs if sources_IDs else []
#         self.parent_consensus = []
#         self.children = []
#
#     def __str__(self):
#         return Sequence.__str__(self) + """ compatibility_to_sources: {0}""".format(  self.compatibility_to_sources)
#
#     def __eq__(self, other):
#         return Sequence.__eq__(self, other) and self.compatibility_to_sources == other.compatibility_to_sources