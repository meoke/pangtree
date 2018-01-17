# class POAGraphRef(object):
#     def __init__(self, sourcesIDs=None, consensus=None):
#         self.sourcesIDs = sourcesIDs if sourcesIDs else []
#         self.consensus = consensus
#
#     def __eq__(self, other):
#         return (self.sourcesIDs == other.sourcesIDs
#                 and self.consensus == other.consensus)
#
#     def __str__(self):
#         return """  SourcesIDs:\n{0},
#                     Consensus:\n{1}""".format("\n".join([str(srcID) for srcID in self.sourcesIDs]),
#                                          str(self.consensus))
