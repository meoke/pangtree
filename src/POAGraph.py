import time
class POAGraph(object):
    def __init__(self, name, title, version, path, sources = None, consensuses = None, nodes = None):
        self.name = name
        self.title = title
        self.version = version
        self.path = path
        self.sources = sources if sources else []
        self.consensuses = consensuses if consensuses else []
        self.nodes = nodes if nodes else []

    def __eq__(self, other):
        return (self.name == other.name
                and self.title == other.title
                and self.version == other.version
                #and self.path == other.path
                and self.nodes == other.nodes
                and self.sources == other.sources
                and self.consensuses == other.consensuses)

    def __str__(self):
        return """  Name: {0}, Title: {1}, Version: {2}, Path: {3},
                    Sources:\n{4},
                    Consensuses:\n{5},
                    Nodes:\n{6}""".format(self.name,
                                         self.title,
                                         self.version,
                                         self.path,
                                         "\n".join([str(source) for source in self.sources]),
                                         "\n".join([str(consensus) for consensus in self.consensuses]),
                                         "\n".join([str(node) for node in self.nodes]))

    def add_node(self, node):
        if node.currentID >= len(self.nodes):
            self.nodes.append(node)
        else:
            self.nodes[node.currentID] = node

    def add_source(self, source):
        new_source_ID = len(self.sources)
        self.sources.append(source)
        for i, node in enumerate(self.nodes):
            if i in source.nodes_IDs:
                self.nodes[i].sources.add(new_source_ID)

    def add_consensus(self, consensus):
        self.consensuses.append(consensus)
        for node in self.nodes:
            node.consensuses_count += 1

    def generate_po(self):
        def generate_introductory_data(active_sources_IDs, active_nodes_IDs):
            nodes_count = len(active_nodes_IDs)
            sources_count = len(active_sources_IDs)

            return ['VERSION=' + self.version,
                    'NAME=' + self.name,
                    'TITLE=' + self.title,
                    'LENGTH=' + str(nodes_count),
                    'SOURCECOUNT=' + str(sources_count)]

        def generate_source_sequences_data(active_sources_IDs):
            def get_source_info(source):
                return '\n'.join(['SOURCENAME=' + source.name,
                                  ' '.join(['SOURCEINFO=', str(len(source.nodes_IDs)),
                                            str(min(source.nodes_IDs)),
                                            str(source.weight),
                                            str(source.consensusID),
                                            str(source.title)])
                                  ])
            self._calc_sources_weights()
            return [get_source_info(self.sources[src_ID]) for src_ID in active_sources_IDs]

        def generate_nodes_data(active_nodes_IDs):
            def get_aligned_nodes_info(node):
                sorted_aligned_nodes = sorted([self.nodes[aligned_node_ID].currentID for aligned_node_ID in node.aligned_to if self.nodes[aligned_node_ID].currentID != -1])
                if sorted_aligned_nodes:
                    if node.currentID > sorted_aligned_nodes[-1]:
                        to_return = "A" + str(sorted_aligned_nodes[0])
                        return to_return
                    to_return = "A" + str(next(node_id for node_id in sorted_aligned_nodes if node_id > node.currentID))
                else:
                    to_return =""
                return to_return
            #
            # def get_sources_info(node_ID):
            #     r = [self.sources[src_ID].currentID for src_ID in self.nodes[node_ID].sources if self.sources[src_ID].currentID is not -1]
            #     return r

            def get_node_info(i, node, nodes_count):
                print("\r\t\tNode " + str(i + 1) + '/' + str(nodes_count), end='')
                l_to_return = ['L' + str(self.nodes[in_node_ID].currentID) for in_node_ID in node.in_nodes if in_node_ID in active_nodes_IDs]
                to_return = "".join([node.base, ":",
                                "".join(l_to_return),
                                "".join(['S' + str(self.sources[src_ID].currentID) for src_ID in node.sources if self.sources[src_ID].active]),
                                get_aligned_nodes_info(node)])

                return to_return

            return [get_node_info(i, node, len(self.nodes)) for i, node in enumerate(self.nodes) if node.currentID != -1]

        active_sources_IDs = [i for i, src in enumerate(self.sources) if src.active]
        active_nodes_IDs = [i for i, node in enumerate(self.nodes) if node.currentID != -1]
        #by default no consensuses and any connected data is printed
        po_lines = []

        po_lines += (generate_introductory_data(active_sources_IDs, active_nodes_IDs))
        po_lines += (generate_source_sequences_data(active_sources_IDs))
        po_lines += (generate_nodes_data(active_nodes_IDs))

        return '\n'.join(po_lines)


    def calculate_compatibility_to_consensuses(self):
        def get_compatibility(source, consensus):
            common_nodes_count = len(set(source.nodes_IDs) & set(consensus.nodes_IDs))
            source_nodes_count = len(source.nodes_IDs)
            return round(common_nodes_count/source_nodes_count,4)

        for consensus in self.consensuses:
            consensus.compatibility_to_sources = [get_compatibility(source, consensus) for source in self.sources]

    def _calc_sources_weights(self):
        def mean(numbers):
            return float(sum(numbers)) / max(len(numbers), 1)

        def get_source_weight(source):
            if not source.active:
                return -1
            else:
                return mean([len(self.nodes[node_ID].sources) for node_ID in source.nodes_IDs])

        def normalize_weight(weight, max_weight, min_weight):
            if weight == -1:
                return -1
            if max_weight - min_weight == 0:
                return 1
            return int(float(format(round((weight - min_weight) / (max_weight - min_weight), 2), '.2f')) * 100)

        weights = [*map(lambda source: get_source_weight(source), self.sources)]
        max_weight = max(weights)
        min_weight = min(set(weights) - set([-1]))
        normalized_weights = [*map(lambda weight: normalize_weight(weight, max_weight, min_weight), weights)]

        for i, source in enumerate(self.sources):
            source.weight = normalized_weights[i]

    def deactivate_different_then(self, sources_currentIDs_to_stay_activate):
        deactivated_sources = 0
        source_current_ID_to_global_ID = {}
        for source_global_ID, source in enumerate(self.sources):
            if source.currentID in sources_currentIDs_to_stay_activate:
                self.sources[source_global_ID].active = True
                self.sources[source_global_ID].currentID = source_global_ID - deactivated_sources
            else:
                self.sources[source_global_ID].active = False
                self.sources[source_global_ID].currentID = -1
                deactivated_sources += 1
                for node_ID in source.nodes_IDs:
                    try:
                        self.nodes[node_ID].sources.remove(source_global_ID)
                    except:
                        pass
            source_current_ID_to_global_ID[self.sources[source_global_ID].currentID] = source_global_ID

        deactivated_nodes = 0
        node_currentID_to_global_ID = {}
        current_node_ID = 0
        for node_global_ID, node in enumerate(self.nodes):
            if len(node.sources) == 0:
                self.nodes[node_global_ID].currentID = -1
                deactivated_nodes += 1
            else:
                #self.nodes[node_global_ID].currentID = node.currentID - deactivated_nodes
                self.nodes[node_global_ID].currentID = current_node_ID
                current_node_ID += 1
            node_currentID_to_global_ID[self.nodes[node_global_ID].currentID] = node_global_ID
        return (source_current_ID_to_global_ID, node_currentID_to_global_ID)


    def activate_sources_with_consensus_unassigned(self):
        currentID=0
        for source_globalID, source in enumerate(self.sources):
            if source.consensusID == -1:
                self.sources[source_globalID].active = True
                self.sources[source_globalID].currentID = currentID
                currentID += 1
                for node_globalID in source.nodes_IDs:
                    self.nodes[node_globalID].sources.add(source_globalID)
            else:
                self.sources[source_globalID].active = False
                self.sources[source_globalID].currentID = -1
                for node_globalID in source.nodes_IDs:
                    try:
                        self.nodes[node_globalID].sources.remove(source_globalID)
                    except:
                        pass

        active_nodes_count = 0
        for node_global_ID, node in enumerate(self.nodes):
            if len(node.sources) == 0:
                self.nodes[node_global_ID].currentID = -1
            else:
                self.nodes[node_global_ID].currentID = active_nodes_count
                active_nodes_count += 1



