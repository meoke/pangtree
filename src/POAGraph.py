class POAGraph(object):
    def __init__(self, name, title, version, path, sources = [], consensuses = [], nodes = []):
        self.name = name
        self.title = title
        self.version = version,
        self.path = path,
        self.sources = sources
        self.consensuses = consensuses
        self.nodes = nodes


    def __eq__(self, other):
        return (self.name == other.name
                and self.title == other.title
                and self.version == other.version
                #and self.path == other.path
                and self.nodedict == other.nodedict
                and self.sources == other.sources
                and self.consensuses == other.consensuses)

    def add_node(self, node):
        self.nodes.append(node)

    def add_source(self, source):
        self.sources.append(source)
