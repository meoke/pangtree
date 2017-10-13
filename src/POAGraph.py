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
