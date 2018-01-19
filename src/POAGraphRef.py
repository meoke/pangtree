import numpy as np

class POAGraphRef(object):
    def __init__(self, sources_IDs=np.array([]), consensus=None, min_compatibility = 0):
        self.sources_IDs = sources_IDs if sources_IDs.size else np.zeros(([]),dtype=np.uint32)
        self.consensus = consensus
        self.min_compatibility = min_compatibility

    def __eq__(self, other):
        return (np.array_equal(self.sourcesIDs, other.sourcesIDs)
                and self.consensus == other.consensus
                and self.min_compatibility == other.min_compatibility)

    def __str__(self):
        return """  SourcesIDs:\n{0},
                    Consensus:\n{1},
                    Min. Compatibility""".format("\n".join([str(srcID) for srcID in self.sourcesIDs]),
                                                            str(self.consensus),
                                                            str(self.min_compatbility))
