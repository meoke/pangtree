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
        if node.ID >= len(self.nodes):
            self.nodes.append(node)
        else:
            self.nodes[node.ID] = node


    def add_source(self, source):
        self.sources.append(source)


    def generate_po(self):
        def generate_introductory_data():
            nodes_count = len([*filter(lambda node : node.sources_count > 0, self.nodes)])
            sources_count = len([*filter(lambda source: source.active, self.sources)])

            return ['VERSION=' + self.version,
                    'NAME=' + self.name,
                    'TITLE=' + self.title,
                    'LENGTH=' + str(nodes_count),
                    'SOURCECOUNT=' + str(sources_count)]

        def generate_source_sequences_data():
            def get_source_info(source):
                return '\n'.join(['SOURCENAME=' + source.name,
                                  ' '.join(['SOURCEINFO=', str(len(source.nodes_IDs)),
                                            str(min(source.nodes_IDs)),
                                            str(source.weight),
                                            str(source.consensusID),
                                            str(source.title)])
                                  ])

            self._calc_sources_weights()
            return [*map(lambda active_source: get_source_info(active_source),
                         filter(lambda source: source.active, self.sources))]

        def generate_nodes_data():
            def get_aligned_nodes_info(node):
                sorted_aligned_nodes = sorted([self.nodes[aligned_node_ID].ID for aligned_node_ID in node.aligned_to if self.nodes[aligned_node_ID].ID != -1])
                if sorted_aligned_nodes:
                    if node.ID > sorted_aligned_nodes[-1]:
                        return "A" + str(sorted_aligned_nodes[0])
                    return "A" + str(next(node_id for node_id in sorted_aligned_nodes if node_id > node.ID))
                else:
                    return ""


            def get_sources_info(node_ID):
                return [i for i, source in enumerate(self.sources) if node_ID in source.nodes_IDs]

            def get_node_info(i, node):
                return "".join([node.base, ":",
                                "".join(['L' + str(self.nodes[in_node_ID].ID) for in_node_ID in node.in_nodes if self.nodes[in_node_ID].ID != -1]),
                                "".join(['S' + str(self.sources[src_ID].ID) for src_ID in get_sources_info(i)]),
                                get_aligned_nodes_info(node)])

            #self._calc_nodes_IDs()

            return [get_node_info(i, node) for i, node in enumerate(self.nodes) if node.ID != -1]


        #by default no consensuses and any connected data is printed
        po_lines = []

        po_lines += (generate_introductory_data())
        po_lines += (generate_source_sequences_data())
        po_lines += (generate_nodes_data())

        return '\n'.join(po_lines)


    def _calculate_compatibility_to_consensuses(self):
        raise Exception("Not implemented")


    def _calc_sources_weights(self):
        def mean(numbers):
            return float(sum(numbers)) / max(len(numbers), 1)

        def get_source_weight(source):
            if not source.active:
                return -1
            else:
                return mean([self.nodes[node_ID].sources_count for node_ID in source.nodes_IDs])

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


    def _calc_nodes_IDs(self):
        raise Exception("Not implemented")
