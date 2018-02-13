import json
import toolkit as t

WEBPAGE_DIR = t.get_real_path('../visualization_template')

class JSONNode:
    def __init__(self, ID, base, source, weight, pos_x, pos_y):
        self.ID = ID
        self.base = base
        self.source = source
        self.weight = weight
        self.pos_x = pos_x
        self.pos_y = pos_y


class JSONEdge:
    def __init__(self, ID, src, target, weight, classes):
        self.ID = ID
        self.src = src
        self.target = target
        self.weight = weight
        self.classes = classes


class JSONConsensus:
    def __init__(self, ID, name, title, sources_comp, length, sources):
        self.ID = ID
        self.name = name
        self.title = title
        self.sources_comp = sources_comp
        self.length = length
        self.sources = sources


class JSONSource:
    def __init__(self, ID, name, title, group, bundles, weight, first_node_ID, length):
        self.ID = ID
        self.name = name
        self.title = title
        self.group = group
        self.bundles = bundles
        self.weight = weight
        self.first_node_ID = first_node_ID
        self.length = length


def _create_common_files(multialignment, output_dir, processing_time):
    # todo skąd się bierze dodatkowy folder o nazwie multialignmentu?
    # wspolne pliki i informacje o multialignmencie (tu będą też bloki? grafy blokowe?)
    common_files_destination = t.join_path(output_dir, 'assets')
    common_files_source = t.join_path(WEBPAGE_DIR, 'assets')
    t.copy_dir(common_files_source, common_files_destination)

    index_template_path = t.join_path(WEBPAGE_DIR, 'index.html')
    with open(index_template_path) as template_file:
        index_content = template_file.read()
        index_content = index_content.replace('info.js', str(multialignment.name) + """_info.js""")
        index_path = t.join_path(output_dir, 'index.html')

    with open(index_path, 'w') as output:
        output.write(index_content)
    #
    #     index_template_path = t.join_path(self.webpage_dir, 'index.html')
    #     with open(index_template_path) as template_file:
    #         index_content = template_file.read()
    #         if consensuses_comparison:
    #             index_content = index_content.replace('sources.js', str(self.poagraph.name) + """_sources_data.js""")
    #             index_content = index_content.replace('consensuses.js', str(self.poagraph.name) + """_consensuses_data.js""")
    #         if graph_visualization:
    #             index_content = index_content.replace('poagraph.js', str(self.poagraph.name) + """_poagraph_data.js""")
    #         index_content = index_content.replace('info.js', str(self.poagraph.name) + """_info.js""")
    #
    #     index_path = t.join_path(self.output_dir, 'index.html')
    #     with open(index_path, 'w') as output:
    #         output.write(index_content)


def _create_poagraph_sources_files(p, output_dir):
    # zawsze - lista sources
    pass


def _create_poagraph_consensus_files(poagraph, output_dir):
    # jeśli consensus option - informacje o consensusach
    pass


def _create_poagraph_graph_files(poagraph, output_dir):
    # jeśli draw_poagraph_option - informacje o poagrafie
    pass


def generate_visualization(multialignment, output_dir, consensus_option, draw_poagraph_option, processing_time):
    _create_common_files(multialignment, output_dir, processing_time)

    for p in multialignment.poagraphs:
        _create_poagraph_sources_files(p, output_dir)

        if consensus_option:
            _create_poagraph_consensus_files(p, output_dir)

        if draw_poagraph_option:
            _create_poagraph_graph_files(p, output_dir)

