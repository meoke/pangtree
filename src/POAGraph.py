import numpy as np
from POAGraphRef import POAGraphRef


class POAGraph(object):
    def __init__(self, name, title, version, path, data_type=-1, sources=None, consensuses=None, nodes=None, ns=np.array([]), nc=np.array([])):
        self.name = name
        self.title = title
        self.version = version
        self.path = path
        self.data_type = data_type
        self.sources = sources if sources else []
        self.consensuses = consensuses if consensuses else []
        self.nodes = nodes if nodes else []
        self.ns = ns if ns.size else np.array([])
        self.nc = nc if nc.size else np.array([])
        self._poagraphrefs = []

    def __eq__(self, other):
        return (self.name == other.name
                and self.title == other.title
                and self.version == other.version
                and self.nodes == other.nodes
                and self.sources == other.sources
                and self.consensuses == other.consensuses
                and np.array_equal(self.ns, other.ns)
                and np.array_equal(self.ns, other.ns))

    def __str__(self):
        return """  Name: {0}, Title: {1}, Version: {2}, Path: {3},
                    Sources:\n{4},
                    Consensuses:\n{5},
                    Nodes:\n{6},
                    NS:\n{7},
                    NC:\n{8}""".format(self.name,
                                         self.title,
                                         self.version,
                                         self.path,
                                         "\n".join([str(source) for source in self.sources]),
                                         "\n".join([str(consensus) for consensus in self.consensuses]),
                                         "\n".join([str(node) for node in self.nodes]),
                                         str(self.ns),
                                         str(self.nc))

    def add_source(self, source): #todo dołożyć jakąś walidację
        self.sources.append(source)

    def add_poagraphref(self, poagrahref, poagraph_parentID):
        poagrahref.ID = len(self._poagraphrefs)
        self._poagraphrefs.append(poagrahref)
        self._poagraphrefs[poagraph_parentID].children_IDs.append(poagrahref.ID)
        return poagrahref.ID

    def add_consensus(self, consensus, consensus_nodes):
        consensus.ID = len(self.consensuses)
        consensus.compatibility_to_sources = np.array([self.get_comp(consensus_nodes, src_ID) for src_ID in range(len(self.sources))])
        self.consensuses.append(consensus)
        if not self.nc.shape[0]:
            self.nc = np.zeros(shape=(1,len(self.nodes)), dtype=np.bool)
            self.nc[consensus.ID] = consensus_nodes
        else:
            self.nc = np.append(self.nc, [consensus_nodes], axis=0)

    def get_min_cutoff(self, poagraphrefs_IDs):
        poagrahrefs_to_check = [p for p in self._poagraphrefs if p.ID in poagraphrefs_IDs]
        return min([p.min_compatibility for p in poagrahrefs_to_check])

    def get_comp(self, consensus_nodes, source_ID):#:calc_compatibility(consensus, tree_node)
        common_nodes_count = np.sum(self.ns[source_ID][:] & consensus_nodes)
        source_nodes_count = np.sum(self.ns[source_ID][:])
        return round(common_nodes_count / source_nodes_count, 4)

    def get_poagraphref_sources_IDs(self, tree_node_ID):
        return self._poagraphrefs[tree_node_ID].sources_IDs

    def create_root_poagraphref(self):
        root_tree_node = POAGraphRef(parent_ID=-1,sources_IDs=np.array(range(len(self.sources))))
        root_tree_node.ID = 0
        self._poagraphrefs.append(root_tree_node)
        return root_tree_node.ID

    def get_poagraphref_compatibility(self, tree_node_ID):
        return self._poagraphrefs[tree_node_ID].min_compatibility

    def set_sources_weights(self):
        def get_source_weight(source_ID):
            source_nodes = np.nonzero(self.ns[source_ID])[0]
            nodes_sizes = [np.sum(self.ns[:, node_ID]) for node_ID in source_nodes]
            return np.mean(nodes_sizes)

        def normalize_weight(weight, max_weight, min_weight):
            if max_weight - min_weight == 0:
                return 100
            return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)

        weights = [*map(lambda source: get_source_weight(source), range(len(self.sources)))]
        max_weight = max(weights)
        min_weight = min(set(weights) - set([-1]))

        for i, src in enumerate(self.sources):
            src.weight = normalize_weight(weights[i], max_weight, min_weight)