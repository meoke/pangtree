import numpy as np


class Sequence(object):
    def __init__(self, ID, name, title):
        self.ID = ID
        self.name = name
        self.title = title

    def __str__(self):
        return """ID: {0},\t name: {1},\t title: {2}""".format(
                self.ID,
                self.name,
                self.title)

    def __eq__(self, other):
        return (self.ID == other.ID
                and self.name == other.name
                and self.title == other.title)


class Source(Sequence):
    def __init__(self, ID, name, title, weight = -1, consensus_ID = -1):
        Sequence.__init__(self, ID=ID, name=name, title=title)
        self.weight = weight
        self.consensus_ID = consensus_ID

    def __str__(self):
        return Sequence.__str__(self) + """\tweight: {0},\tconsensus_ID: {1}""".format(
            self.weight, self.consensus_ID)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and \
            self.weight == other.weight and \
            self.consensus_ID == other.consensus_ID


class Consensus(Sequence):
    def __init__(self, ID, name, title, compatibility_to_sources=np.array([])):
        Sequence.__init__(self, ID=ID, name=name, title=title)
        self.compatibility_to_sources = compatibility_to_sources if compatibility_to_sources.size else np.array([], dtype=int)

    def __str__(self):
        return Sequence.__str__(self) + """ compatibility_to_sources: {0}""".format(  self.compatibility_to_sources)

    def __eq__(self, other):
        return Sequence.__eq__(self, other) and np.array_equal(self.compatibility_to_sources, other.compatibility_to_sources)