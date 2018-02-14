import numpy as np
import toolkit as t
from POAGraphRef import POAGraphRef
# import consensus as cons

class POAGraph(object):
    def __init__(self, name, title, version, path, data_type=-1, sources=None, consensuses=None, nodes=None, ns=np.array([]), nc=np.array([])):
        self.name = name
        self.title = title
        self.version = version
        self.path = path #todo czy potrzebne? A może zarządzane przez Multialignment?
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
                #and self.path == other.path
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

    # def clean(self):
    #     for source in self.sources:
    #         source.truncate_nodes_ids()

    # def add_node(self, node):
    #     self.nodes.append(node)
        # if node.currentID >= len(self.nodes):
        #     self.nodes.append(node)
        # else:
        #     self.nodes[node.currentID] = node

    def add_source(self, source): #todo dołożyć jakąś walidację
        # new_source_ID = len(self.sources)
        self.sources.append(source)
        # for i, node in enumerate(self.nodes):
        #     if i in source.nodes_IDs:
        #         self.nodes[i].sources.add(new_source_ID)

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
        # np.reshape(self.nc, newshape=(consensus.ID +1, len(self.nodes)))
        # self.nc[consensus.ID] = consensus_nodes

    def get_min_cutoff(self, poagraphrefs_IDs):
        poagrahrefs_to_check = [p for p in self._poagraphrefs if p.ID in poagraphrefs_IDs]
        return min([p.min_compatibility for p in poagrahrefs_to_check])

    def get_comp(self, consensus_nodes, source_ID):#:calc_compatibility(consensus, tree_node)
        common_nodes_count = sum(self.ns[source_ID][:] & consensus_nodes)
        source_nodes_count = sum(self.ns[source_ID][:])
        return round(common_nodes_count / source_nodes_count, 4)

        # def get_compatibility(source, consensus):
        #     common_nodes_count = len(set(source.nodes_IDs) & set(consensus.nodes_IDs))
        #     source_nodes_count = len(source.nodes_IDs)
        #     return round(common_nodes_count / source_nodes_count, 4)

        # if not consensusID:
        #     consensuses_to_calculate = [consensus for consensus in self.consensuses if consensus.level == level]
        # else:
        #     consensuses_to_calculate = [self.consensuses[consensusID]]

        # if not sources_IDs:
        #     srcs = [source for source in self.sources]
        # else:
        #     srcs = [src for src in self.sources if src.currentID in sources_IDs]

        # for consensus in consensuses_to_calculate:
        #     # consensus.compatibility_to_sources = [get_compatibility(source, consensus) for source in self.sources]
        #     consensus.compatibility_to_sources = [get_compatibility(source, consensus) for source in srcs]


    # teraz to robi po_writer
    # def save_as_po(self, output_dir, sources_IDs):
    #     def write_introduction_data(output_po_file, active_nodes_count):
    #         output_po_file.writelines('VERSION=' + self.version + '\n')
    #         output_po_file.writelines('NAME=' + self.name + '\n')
    #         output_po_file.writelines('TITLE=' + self.title + '\n')
    #         output_po_file.writelines('LENGTH=' + str(active_nodes_count) + '\n')
    #         output_po_file.writelines('SOURCECOUNT=' + str(len(sources_IDs)) + '\n')
    #
    #     def write_source_sequences(output_po_file, nodes):
    #         def get_source_info(source):
    #             return ['SOURCENAME=' + source.name,
    #                     '\n',
    #                     " ".join(['SOURCEINFO=', str(len(source.nodes_IDs)),
    #                     str((nodes['temp_nodeID'][nodes['org_nodeID'] == source.nodes_IDs[0]])[0]),
    #                     str(source.weight),
    #                     str(-1),
    #                     str(source.title)]),
    #                     '\n']
    #
    #         self._calc_partial_sources_weights(sources_IDs, nodes)
    #         for srcID in sources_IDs:
    #             output_po_file.writelines(get_source_info(self.sources[srcID]))
    #
    #     def write_nodes(output_po_file, nodes, sources):
    #         def get_node_info(node):
    #             L_to_return = ['L' + str((nodes['temp_nodeID'][nodes['org_nodeID'] == in_nodeID])[0]) for in_nodeID in node.in_nodes
    #                            if in_nodeID in nodes['org_nodeID'][nodes['active']==True]]
    #             S_to_return = ['S' + str((sources['temp_srcID'][sources['org_srcID']== srcID])[0]) for srcID in sources_IDs
    #                            if node.ID in self.sources[srcID].nodes_IDs]
    #             if node.aligned_to:
    #                 if node.aligned_to in nodes['org_nodeID'][nodes['active']==True]:
    #                     A_to_return = "A" + str(nodes['temp_nodeID'][nodes['org_nodeID'] == node.aligned_to][0])
    #                 else:
    #                     A_to_return = ""
    #             else:
    #                 A_to_return = ""
    #
    #             return [node.base, ':',
    #                     "".join(L_to_return),
    #                     "".join(S_to_return),
    #                     A_to_return,
    #                     "\n"
    #                     ]
    #
    #
    #         for node in self.nodes:
    #             if node.ID in np.where(nodes['active'] == True)[0]:
    #                 output_po_file.writelines(get_node_info(node))
    #
    #     nodes = np.zeros(shape=(len(self.nodes)), dtype = [('org_nodeID', np.uint32),
    #                                                        ('temp_nodeID', np.uint32),
    #                                                        ('active', np.bool),
    #                                                        ('sources_count', np.uint16)])
    #     for i, src in enumerate(self.sources):
    #         nodes['active'][src.nodes_IDs] = True
    #         nodes['sources_count'][src.nodes_IDs] = nodes['sources_count'][src.nodes_IDs] + 1
    #
    #     active_nodes_count = len(nodes[nodes['active'] == True])
    #     nodes['org_nodeID'][nodes['active'] == True] = np.where(nodes['active'] == True)[0]
    #     nodes['temp_nodeID'][nodes['active'] == True] = range(len(nodes[nodes['active']==True]))
    #
    #
    #
    #     sources = np.zeros(shape=(len(self.sources)), dtype = [('org_srcID', np.uint32),
    #                                                            ('temp_srcID', np.uint32),
    #                                                            ('active', np.bool)])
    #
    #     sources['active'][sources_IDs] = True
    #     sources['org_srcID'][sources['active'] == True] = np.where(sources['active'] == True)[0]
    #     sources['temp_srcID'][sources['active'] == True] = range(len(sources[sources['active']==True]))
    #
    #     po_file_name = t.get_next_child_file_name(output_dir, self.name + str(".po"))
    #     with open(po_file_name, 'w') as output_po_file:
    #         write_introduction_data(output_po_file, active_nodes_count)
    #         write_source_sequences(output_po_file, nodes)
    #         write_nodes(output_po_file, nodes, sources)
    #
    #     return po_file_name, nodes
    # output_po_file.write(poagraph_as_po)

    def get_poagraphref_sources_IDs(self, tree_node_ID):
        return self._poagraphrefs[tree_node_ID].sources_IDs

    def create_parent_poagraphref(self):
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
            # if weight == -1:
            #     return -1
            if max_weight - min_weight == 0:
                return 100
            return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)

        weights = [*map(lambda source: get_source_weight(source), range(len(self.sources)))]
        max_weight = max(weights)
        min_weight = min(set(weights) - set([-1]))

        for i, src in enumerate(self.sources):
            src.weight = normalize_weight(weights[i], max_weight, min_weight)

#     def generate_source_sequences_data(active_sources_IDs):
#         def get_source_info(source):
#             return '\n'.join(['SOURCENAME=' + source.name,
#                               ' '.join(['SOURCEINFO=', str(len(source.nodes_IDs)),
#                                         str(min(source.nodes_IDs)),
#                                         str(source.weight),
#                                         str(source.consensusID),
#                                         str(source.title)])
#                               ])
#         self._calc_sources_weights()
#         return [get_source_info(self.sources[src_ID]) for src_ID in active_sources_IDs]
#
#     def generate_nodes_data(active_nodes_IDs):
#         def get_aligned_nodes_info(node):
#             sorted_aligned_nodes = sorted([self.nodes[aligned_node_ID].currentID for aligned_node_ID in node.aligned_to if self.nodes[aligned_node_ID].currentID != -1])
#             if sorted_aligned_nodes:
#                 if node.currentID > sorted_aligned_nodes[-1]:
#                     to_return = "A" + str(sorted_aligned_nodes[0])
#                     return to_return
#                 to_return = "A" + str(next(node_id for node_id in sorted_aligned_nodes if node_id > node.currentID))
#             else:
#                 to_return =""
#             return to_return
#
#         def get_node_info(i, node, nodes_count):
#             print("\r\t\tNode " + str(i + 1) + '/' + str(nodes_count), end='')
#             l_to_return = ['L' + str(self.nodes[in_node_ID].currentID) for in_node_ID in node.in_nodes if in_node_ID in active_nodes_IDs]
#             to_return = "".join([node.base, ":",
#                             "".join(l_to_return),
#                             "".join(['S' + str(self.sources[src_ID].currentID) for src_ID in node.sources if self.sources[src_ID].active]),
#                             get_aligned_nodes_info(node)])
#
#             return to_return
#
#         return [get_node_info(i, node, len(self.nodes)) for i, node in enumerate(self.nodes) if node.currentID != -1]
#
#     active_sources_IDs = [i for i, src in enumerate(self.sources) if src.active]
#     active_nodes_IDs = [i for i, node in enumerate(self.nodes) if node.currentID != -1]
#     #by default no consensuses and any connected data is printed
#     po_lines = []
#
#     po_lines += (generate_introductory_data(active_sources_IDs, active_nodes_IDs))
#     po_lines += (generate_source_sequences_data(active_sources_IDs))
#     po_lines += (generate_nodes_data(active_nodes_IDs))
#
#     return '\n'.join(po_lines)

    # def get_partial_sources_weights(self, sources_IDs, nodes):
    #     def get_source_weight(source):
    #         # if not source.ID in sources_IDs:
    #         #     return -1
    #         # else:
    #         src_nodes_IDs = self.ns[source] == True
    #         # return np.mean(np.array([nodes['sources_count'][nodes['org_nodeID'] == node_ID][0] for node_ID in source.nodes_IDs]))
    #         return np.mean(np.array([nodes['sources_count'][nodes['orig_ID'] == node_ID][0] for node_ID in src_nodes_IDs]))
    #
    #     def normalize_weight(weight, max_weight, min_weight):
    #         if weight == -1:
    #             return -1
    #         if max_weight - min_weight == 0:
    #             return 1
    #         return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)
    #
    #     weights = [*map(lambda source: get_source_weight(source), sources_IDs)]
    #     max_weight = max(weights)
    #     min_weight = min(set(weights) - set([-1]))
    #
    #     normalized_weights = [*map(lambda weight: normalize_weight(weight, max_weight, min_weight), weights)]
    #     return normalized_weights
        # for i, source in enumerate(self.sources):
        #     source.weight = normalized_weights[i]

    # zakomentowane 6.02 20:20
    # def _calc_partial_sources_weights(self, sourcesIDs, nodes):
    #     def get_source_weight(source):
    #         if not source.ID in sourcesIDs:
    #             return -1
    #         else:
    #             return np.mean(np.array([nodes['sources_count'][nodes['org_nodeID'] == node_ID][0] for node_ID in source.nodes_IDs]))
    #
    #     def normalize_weight(weight, max_weight, min_weight):
    #         if weight == -1:
    #             return -1
    #         if max_weight - min_weight == 0:
    #             return 1
    #         return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)
    #
    #     weights = [*map(lambda source: get_source_weight(source), self.sources)]
    #     max_weight = max(weights)
    #     min_weight = min(set(weights) - set([-1]))
    #
    #     normalized_weights = [*map(lambda weight: normalize_weight(weight, max_weight, min_weight), weights)]
    #     for i, source in enumerate(self.sources):
    #         source.weight = normalized_weights[i]

        # def _calc_partial_sources_weights(self, sourcesIDs_to_use, new_to_original_nodes_IDs):
        #     def mean(numbers):
        #         return float(sum(numbers)) / max(len(numbers), 1)
        #
        #     def get_source_weight(source, nodes_sources_count):
        #         if not source.currentID in sourcesIDs_to_use:
        #             return -1
        #         else:
        #             return mean([nodes_sources_count[node_ID] for node_ID in source.nodes_IDs])
        #
        #     def normalize_weight(weight, max_weight, min_weight):
        #         if weight == -1:
        #             return -1
        #         if max_weight - min_weight == 0:
        #             return 1
        #         return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)
        #
        #
        #     nodes_sources_count = {}
        #     for new_nodeID, org_nodeID in new_to_original_nodes_IDs.items():
        #         nodes_sources_count[org_nodeID] = len([srcID for srcID in self.nodes[org_nodeID].sources if srcID in sourcesIDs_to_use])
        #
        #     weights = [*map(lambda source: get_source_weight(source, nodes_sources_count), self.sources)]
        #     max_weight = max(weights)
        #     min_weight = min(set(weights) - set([-1]))
        #
        #     normalized_weights = [*map(lambda weight: normalize_weight(weight, max_weight, min_weight), weights)]
        #
        #     for i, source in enumerate(self.sources):
        #         source.weight = normalized_weights[i]
        #

    # def run_tree_consensus_generation(self):
    #     consensus_output_dir = t.create_child_dir(self.path, "tconsensus")
    #     hbmin = 0.2
    #
    #     parent_tree_node = POAGraphRef(sources_IDs=np.array(range(len(self.sources))))
    #     tree_nodes_to_process = [parent_tree_node]
    #
    #     finished_sources_IDs = np.empty(len(self.sources), dtype=np.int16)
    #     finished_sources_IDs.fill(-1)
    #     while -1 in finished_sources_IDs:
    #         new_tree_nodes_to_process = []
    #         for tree_node in tree_nodes_to_process:
    #             new_tree_nodes_to_process = cons.process_tree_node(tree_node, consensus_output_dir, consensus_output_dir)


    # def add_consensus(self, consensus):
    #     self.consensuses.append(consensus)
    #     #TODO to chyba jest źle i niepotrzebne...
    #     for node in self.nodes:
    #         node.consensuses_count += 1
    #
    # def remove_last_consensus(self):
    #     self.consensuses.pop()
    #     # TODO to chyba jest źle i niepotrzebne...
    #     for node in self.nodes:
    #         node.consensuses_count -= 1
    #
    # def generate_po(self):
    #     def generate_introductory_data(active_sources_IDs, active_nodes_IDs):
    #         nodes_count = len(active_nodes_IDs)
    #         sources_count = len(active_sources_IDs)
    #
    #         return ['VERSION=' + self.version,
    #                 'NAME=' + self.name,
    #                 'TITLE=' + self.title,
    #                 'LENGTH=' + str(nodes_count),
    #                 'SOURCECOUNT=' + str(sources_count)]
    #
    #     def generate_source_sequences_data(active_sources_IDs):
    #         def get_source_info(source):
    #             return '\n'.join(['SOURCENAME=' + source.name,
    #                               ' '.join(['SOURCEINFO=', str(len(source.nodes_IDs)),
    #                                         str(min(source.nodes_IDs)),
    #                                         str(source.weight),
    #                                         str(source.consensusID),
    #                                         str(source.title)])
    #                               ])
    #         self._calc_sources_weights()
    #         return [get_source_info(self.sources[src_ID]) for src_ID in active_sources_IDs]
    #
    #     def generate_nodes_data(active_nodes_IDs):
    #         def get_aligned_nodes_info(node):
    #             sorted_aligned_nodes = sorted([self.nodes[aligned_node_ID].currentID for aligned_node_ID in node.aligned_to if self.nodes[aligned_node_ID].currentID != -1])
    #             if sorted_aligned_nodes:
    #                 if node.currentID > sorted_aligned_nodes[-1]:
    #                     to_return = "A" + str(sorted_aligned_nodes[0])
    #                     return to_return
    #                 to_return = "A" + str(next(node_id for node_id in sorted_aligned_nodes if node_id > node.currentID))
    #             else:
    #                 to_return =""
    #             return to_return
    #
    #         def get_node_info(i, node, nodes_count):
    #             print("\r\t\tNode " + str(i + 1) + '/' + str(nodes_count), end='')
    #             l_to_return = ['L' + str(self.nodes[in_node_ID].currentID) for in_node_ID in node.in_nodes if in_node_ID in active_nodes_IDs]
    #             to_return = "".join([node.base, ":",
    #                             "".join(l_to_return),
    #                             "".join(['S' + str(self.sources[src_ID].currentID) for src_ID in node.sources if self.sources[src_ID].active]),
    #                             get_aligned_nodes_info(node)])
    #
    #             return to_return
    #
    #         return [get_node_info(i, node, len(self.nodes)) for i, node in enumerate(self.nodes) if node.currentID != -1]
    #
    #     active_sources_IDs = [i for i, src in enumerate(self.sources) if src.active]
    #     active_nodes_IDs = [i for i, node in enumerate(self.nodes) if node.currentID != -1]
    #     #by default no consensuses and any connected data is printed
    #     po_lines = []
    #
    #     po_lines += (generate_introductory_data(active_sources_IDs, active_nodes_IDs))
    #     po_lines += (generate_source_sequences_data(active_sources_IDs))
    #     po_lines += (generate_nodes_data(active_nodes_IDs))
    #
    #     return '\n'.join(po_lines)
    #
    # def generate_partial_po(self, sourcesIDs_to_use):
    #     def generate_introductory_data(nodes_IDs):
    #         nodes_count = len(nodes_IDs)
    #         sources_count = len(sourcesIDs_to_use)
    #
    #         return ['VERSION=' + self.version,
    #                 'NAME=' + self.name,
    #                 'TITLE=' + self.title,
    #                 'LENGTH=' + str(nodes_count),
    #                 'SOURCECOUNT=' + str(sources_count)]
    #
    #     def generate_source_sequences_data(new_to_original_nodes_IDs):
    #         def get_source_info(source):
    #             return '\n'.join(['SOURCENAME=' + source.name,
    #                               ' '.join(['SOURCEINFO=', str(len(source.nodes_IDs)),
    #                                         str(min(new_to_original_nodes_IDs.keys())),
    #                                         str(source.weight),
    #                                         str(source.consensusID),
    #                                         str(source.title)])
    #                               ])
    #         self._calc_partial_sources_weights(sourcesIDs_to_use, new_to_original_nodes_IDs)
    #         return [get_source_info(self.sources[src_ID]) for src_ID in sourcesIDs_to_use]
    #
    #     def generate_nodes_data(new_to_original_nodes_IDs, original_to_new_nodes_IDs):
    #         def get_aligned_nodes_info(node):
    #             sorted_aligned_nodes = sorted([original_to_new_nodes_IDs[self.nodes[aligned_node_ID].currentID]
    #                                            for aligned_node_ID in node.aligned_to if self.nodes[aligned_node_ID].currentID in original_to_new_nodes_IDs.keys()])
    #             if sorted_aligned_nodes:
    #                 if original_to_new_nodes_IDs[node.currentID] > sorted_aligned_nodes[-1]:
    #                     to_return = "A" + str(sorted_aligned_nodes[0])
    #                     return to_return
    #                 to_return = "A" + str(next(node_id for node_id in sorted_aligned_nodes if node_id > original_to_new_nodes_IDs[node.currentID]))
    #             else:
    #                 to_return =""
    #             return to_return
    #
    #         def get_node_info(node, nodes_count, original_to_new_sources_IDs):
    #             print("\r\t\tNode " + str(original_to_new_nodes_IDs[node.currentID] + 1) + '/' + str(nodes_count), end='')
    #             l_to_return = ['L' + str(original_to_new_nodes_IDs[in_node_ID]) for in_node_ID in node.in_nodes if in_node_ID in original_to_new_nodes_IDs.keys()]
    #             to_return = "".join([node.base, ":",
    #                             "".join(l_to_return),
    #                             "".join(['S' + str(original_to_new_sources_IDs[src_ID]) for src_ID in node.sources if src_ID in sourcesIDs_to_use]),
    #                             get_aligned_nodes_info(node)])
    #
    #             return to_return
    #
    #         return [get_node_info(self.nodes[org_node_ID], len(original_to_new_nodes_IDs), original_to_new_sources_IDs) for org_node_ID, new_node_ID in sorted(original_to_new_nodes_IDs.items())]
    #
    #     new_to_original_nodes_IDs = {}
    #     original_to_new_nodes_IDs = {}
    #     new_nodes_count = 0
    #     for i, node in enumerate(self.nodes):
    #         if any([srcID in sourcesIDs_to_use for srcID in node.sources]):
    #             new_to_original_nodes_IDs[new_nodes_count] = node.currentID
    #             original_to_new_nodes_IDs[node.currentID] = new_nodes_count
    #             new_nodes_count +=1
    #
    #     new_to_original_sources_IDs = {}
    #     original_to_new_sources_IDs = {}
    #     new_sources_count = 0
    #     for i, source in enumerate(self.sources):
    #         if source.currentID in sourcesIDs_to_use:
    #             new_to_original_sources_IDs[new_sources_count] = source.currentID
    #             original_to_new_sources_IDs[source.currentID] = new_sources_count
    #             new_sources_count += 1
    #
    #     #by default no consensuses and any connected data is printed
    #     po_lines = []
    #
    #     po_lines += (generate_introductory_data(new_to_original_nodes_IDs.keys()))
    #     po_lines += (generate_source_sequences_data(new_to_original_nodes_IDs))
    #     po_lines += (generate_nodes_data(new_to_original_nodes_IDs, original_to_new_nodes_IDs))
    #
    #     return ('\n'.join(po_lines), new_to_original_nodes_IDs)
    #
    # def calculate_compatibility_to_consensuses(self, consensusID=None, sources_IDs = None, level=-1):
    #     def get_compatibility(source, consensus):
    #         common_nodes_count = len(set(source.nodes_IDs) & set(consensus.nodes_IDs))
    #         source_nodes_count = len(source.nodes_IDs)
    #         return round(common_nodes_count/source_nodes_count,4)
    #
    #     if not consensusID:
    #         consensuses_to_calculate = [consensus for consensus in self.consensuses if consensus.level == level]
    #     else:
    #         consensuses_to_calculate = [self.consensuses[consensusID]]
    #
    #     if not sources_IDs:
    #         srcs = [source for source in self.sources]
    #     else:
    #         srcs = [src for src in self.sources if src.currentID in sources_IDs]
    #
    #     for consensus in consensuses_to_calculate:
    #         # consensus.compatibility_to_sources = [get_compatibility(source, consensus) for source in self.sources]
    #         consensus.compatibility_to_sources = [get_compatibility(source, consensus) for source in srcs]
    #

    #
    # def _calc_sources_weights(self):
    #     def mean(numbers):
    #         return float(sum(numbers)) / max(len(numbers), 1)
    #
    #     def get_source_weight(source):
    #         if not source.active:
    #             return -1
    #         else:
    #             return mean([len(self.nodes[node_ID].sources) for node_ID in source.nodes_IDs])
    #
    #     def normalize_weight(weight, max_weight, min_weight):
    #         if weight == -1:
    #             return -1
    #         if max_weight - min_weight == 0:
    #             return 1
    #         return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)
    #
    #     weights = [*map(lambda source: get_source_weight(source), self.sources)]
    #     max_weight = max(weights)
    #     min_weight = min(set(weights) - set([-1]))
    #     normalized_weights = [*map(lambda weight: normalize_weight(weight, max_weight, min_weight), weights)]
    #
    #     for i, source in enumerate(self.sources):
    #         source.weight = normalized_weights[i]
    #
    # def deactivate_different_then(self, sources_currentIDs_to_stay_activate):
    #     deactivated_sources = 0
    #     source_current_ID_to_global_ID = {}
    #     for source_global_ID, source in enumerate(self.sources):
    #         if source.currentID in sources_currentIDs_to_stay_activate:
    #             self.sources[source_global_ID].active = True
    #             self.sources[source_global_ID].currentID = source_global_ID - deactivated_sources
    #         else:
    #             self.sources[source_global_ID].active = False
    #             self.sources[source_global_ID].currentID = -1
    #             deactivated_sources += 1
    #             for node_ID in source.nodes_IDs:
    #                 try:
    #                     self.nodes[node_ID].sources.remove(source_global_ID)
    #                 except:
    #                     pass
    #         source_current_ID_to_global_ID[self.sources[source_global_ID].currentID] = source_global_ID
    #
    #     deactivated_nodes = 0
    #     node_currentID_to_global_ID = {}
    #     current_node_ID = 0
    #     for node_global_ID, node in enumerate(self.nodes):
    #         if len(node.sources) == 0:
    #             self.nodes[node_global_ID].currentID = -1
    #             deactivated_nodes += 1
    #         else:
    #             self.nodes[node_global_ID].currentID = current_node_ID
    #             current_node_ID += 1
    #         node_currentID_to_global_ID[self.nodes[node_global_ID].currentID] = node_global_ID
    #     return (source_current_ID_to_global_ID, node_currentID_to_global_ID)
    #
    # def activate_sources_with_consensus_unassigned(self):
    #     currentID=0
    #     for source_globalID, source in enumerate(self.sources):
    #         if source.consensusID == -1:
    #             self.sources[source_globalID].active = True
    #             self.sources[source_globalID].currentID = currentID
    #             currentID += 1
    #             for node_globalID in source.nodes_IDs:
    #                 self.nodes[node_globalID].sources.add(source_globalID)
    #         else:
    #             self.sources[source_globalID].active = False
    #             self.sources[source_globalID].currentID = -1
    #             for node_globalID in source.nodes_IDs:
    #                 try:
    #                     self.nodes[node_globalID].sources.remove(source_globalID)
    #                 except:
    #                     pass
    #
    #     active_nodes_count = 0
    #     for node_global_ID, node in enumerate(self.nodes):
    #         if len(node.sources) == 0:
    #             self.nodes[node_global_ID].currentID = -1
    #         else:
    #             self.nodes[node_global_ID].currentID = active_nodes_count
    #             active_nodes_count += 1
    #
    #
    #

