import json
import numpy as np
import toolkit as t
from data_types import ebola as eb
from Sequence import Source

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

class SourceEncoder(json.JSONEncoder):
    poagraph = None

    def get_names_and_group(self, source):
        if self.poagraph.data_type is 'ebola':
            source_name = eb.extract_ebola_code_part(source.name, 0)
            source_alt_name = eb.get_official_ebola_name(eb.extract_ebola_code_part(source.title, 1))
            source_group_name = eb.get_ebola_group_name(eb.extract_ebola_code_part(source.title, 1))
        else:
            source_name = source.name
            source_alt_name = '-'
            source_group_name = '-'
        return source_name, source_alt_name, source_group_name

    def get_first_node_ID(self, source_ID):
        return int(np.nonzero(self.poagraph.ns[source_ID])[0][0])

    def get_length(self, source_ID):
        return int(np.sum(self.poagraph.ns[source_ID]))

    def default(self, obj):
        if isinstance(obj, Source):
            name, alt_name, group_name = self.get_names_and_group(obj)
            bundles = -1  # todo bundles
            weight = obj.weight
            first_node_ID = self.get_first_node_ID(obj.ID)
            length = self.get_length(obj.ID)
            return {"ID": obj.ID,
                    "name": name,
                    "title": alt_name,
                    "group": group_name,
                    "bundle_ID": bundles,
                    "weight": weight,
                    "first_node_ID": first_node_ID,
                    "length": length }
        return json.JSONEncoder.default(self, obj)

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

def _get_poagraph_output_dir(p, output_dir):
    return t.create_next_sibling_dir(output_dir, p.name)

def _create_poagraph_sources_files(p, output_dir):
    def get_names_and_group(source):
        if data_type is 'ebola':
            source_name = eb.extract_ebola_code_part(source.name, 0)
            source_alt_name = eb.get_official_ebola_name(eb.extract_ebola_code_part(source.title, 1))
            source_group_name = eb.get_ebola_group_name(eb.extract_ebola_code_part(source.title, 1))
        else:
            source_name = source.name
            source_alt_name = '-'
            source_group_name = '-'
        return source_name, source_alt_name, source_group_name

    def get_first_node_ID(source_ID):
        return np.nonzero(p.ns[source_ID])[0]

    def get_length(source_ID):
        return np.sum(p.ns[source_ID])

    # JSONSources = []
    # for src in p.sources:
    #     name, alt_name, group_name = get_names_and_group(src)
    #     bundles = -1 #todo bundles
    #     weight = src.weight
    #     first_node_ID = get_first_node_ID(src.ID)
    #     length = get_length(src.ID)
    #     JSONSources.append(JSONSource(src.ID,
    #                                   name,
    #                                   alt_name,
    #                                   group_name,
    #                                   bundles,
    #                                   weight,
    #                                   first_node_ID,
    #                                   length))
    sources_filename = t.join_path(output_dir, "sources.js")
    SourceEncoder.poagraph = p
    l = p.sources
    with open(sources_filename, 'w') as out:
        json.dump(l, fp=out, cls=SourceEncoder)


def _create_poagraph_consensus_files(poagraph, output_dir):
    # jeśli consensus option - informacje o consensusach
    pass


def _create_poagraph_graph_files(poagraph, output_dir):
    # jeśli draw_poagraph_option - informacje o poagrafie
    pass


def generate_visualization(multialignment, output_dir, consensus_option, draw_poagraph_option, processing_time):
    _create_common_files(multialignment, output_dir, processing_time)

    for p in multialignment.poagraphs:
        poagraph_output_dir = _get_poagraph_output_dir(p, output_dir)
        _create_poagraph_sources_files(p, poagraph_output_dir)

        if consensus_option:
            _create_poagraph_consensus_files(p, poagraph_output_dir)

        if draw_poagraph_option:
            _create_poagraph_graph_files(p, poagraph_output_dir)

