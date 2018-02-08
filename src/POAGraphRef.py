import numpy as np

class POAGraphRef(object):
    def __init__(self, sources_IDs=np.array([]), consensus_ID=None, min_compatibility = 0):
        self.ID = -1
        self.sources_IDs = sources_IDs if sources_IDs.size else np.zeros(([]),dtype=np.uint32)
        self.consensus_ID = consensus_ID
        self.min_compatibility = min_compatibility #todo czy to jest min_compatibility to swoich sourcesIDs czy ogólnie? gdzie są ogólnie compatibility zapisane?
        self.children_IDs = []
        self.parent_ID = []

    def __eq__(self, other):
        return (np.array_equal(self.sourcesIDs, other.sourcesIDs)
                and self.consensus_ID == other.consensus_ID
                and self.min_compatibility == other.min_compatibility)

    def __str__(self):
        return """  SourcesIDs:\n{0},
                    Consensus:\n{1},
                    Min. Compatibility""".format("\n".join([str(src_ID) for src_ID in self.sources_IDs]),
                                                 str(self.consensus_ID),
                                                 str(self.min_compatibility))
