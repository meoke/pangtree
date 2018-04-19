import json
import numpy as np
import toolkit as t
from data_types import ebola as eb
from data_types import mycoplasma as myc
from Sequence import Source, Consensus
from POAGraphRef import POAGraphRef
from Node import Node
import maf_reader
from nucleotides import _nucleotide_inverse_dictionary

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


class SourceEncoder(json.JSONEncoder):
    poagraph = None

    def get_names_and_group(self, source):
        if self.poagraph.data_type == 'ebola':
            source_name = eb.extract_ebola_code_part(source.name, 0)
            source_alt_name = eb.get_official_ebola_name(eb.extract_ebola_code_part(source.title, 1))
            source_group_name = eb.get_ebola_group_name(eb.extract_ebola_code_part(source.title, 1))
        elif self.poagraph.data_type == 'mycoplasma':
            source_name = source.name
            source_alt_name = myc.get_strain_name(source.name)
            source_group_name = myc.get_group_name(source.name)
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


class POAGraphRefEncoder(json.JSONEncoder):
    poagraph = None

    def get_length(self, consensus_ID):
        return int(np.sum(self.poagraph.nc[consensus_ID]))

    def default(self, obj):
        if isinstance(obj, POAGraphRef):
            if obj.consensus_ID is None:
                return {"ID": obj.ID,
                        "name": "",
                        "title": "",
                        "sources_compatibility": [],
                        "length": 0,
                        "sources": obj.sources_IDs.tolist(),
                        "parent": obj.parent_ID,
                        "level": obj.min_compatibility,
                        "children": obj.children_IDs}
            consensus = self.poagraph.consensuses[obj.consensus_ID]

            return {"ID": obj.ID,
                    "name": consensus.name,
                    "title": consensus.title,
                    "sources_compatibility": consensus.compatibility_to_sources.tolist(),
                    "length": self.get_length(consensus.ID),
                    "sources": obj.sources_IDs.tolist(),
                    "parent": obj.parent_ID,
                    "level": obj.min_compatibility,
                    "children": obj.children_IDs}

        return json.JSONEncoder.default(self, obj)


class POAGraphNodeEncoder(json.JSONEncoder):
    poagraph = None
    def get_source_list(self, node_ID):
        return np.where(self.poagraph.ns[:, node_ID] == True)[0]

    def get_weight(self, node_ID):
        return 1

    def get_pos(self, node_ID):
        n = self.poagraph.nodes[node_ID]
        if n.aligned_to is None:
            y = 0
            x = node_ID * 20
        else:
            aligned_nodes = []
            while True:
                next_aligned_node = n.aligned_to
                if next_aligned_node in aligned_nodes:
                    break
                aligned_nodes.append(next_aligned_node)
                n = self.poagraph.nodes[n.aligned_to]
            sorted_aligned_nodes = sorted(aligned_nodes)
            node_index = sorted_aligned_nodes.index(node_ID)
            if len(aligned_nodes) == 4:
                positions = [-40, -20, 20, 40]
                y = positions[node_index]
                x = 20*(sorted_aligned_nodes[1] + sorted_aligned_nodes[2])/2
            elif len(aligned_nodes) == 3:
                positions = [-30, 0, 30]
                y = positions[node_index]
                x = 20*(sorted_aligned_nodes[1])
            elif len(aligned_nodes) == 2:
                positions = [-30, 30]
                y = positions[node_index]
                x = 20*(sorted_aligned_nodes[0] + sorted_aligned_nodes[1]) / 2

        return x, y

    def default(self, obj):
        if isinstance(obj, Node):
            sources = self.get_source_list(obj.ID)
            weight = self.get_weight(obj.ID)
            x, y = self.get_pos(obj.ID)
            return {
                      "data": {
                        "id": obj.ID,
                        "nucleobase": obj.base,
                        "source": sources,
                        "weight": weight
                      },
                      "position": { "x": x, "y": y }
                    }
        return {}
        return json.JSONEncoder.default(self, obj)


class Edge(object):
    def __init__(self, id, source, target, weight, consensus, level, classes):
        self.id = id
        self.source = source
        self.target = target
        self.weight = weight
        self.consensus = consensus
        self.level = level
        self.classes = classes


class POAGraphEdgeEncoder(json.JSONEncoder):
    poagraph = None
    def default(self, obj):
        if isinstance(obj, Edge):
            return {
                      "data": {
                        "id": obj.id,
                        "source": obj.source,
                        "target": obj.target,
                        "weight": obj.weight,
                        "consensus": obj.consensus,
                        "level": obj.level
                      },
                      "classes": obj.classes
                    }
        return {}
        return json.JSONEncoder.default(self, obj)



def generate_visualization(multialignment, output_dir, consensus_option, draw_poagraph_option, blocks_option, processing_time):
    _create_common_files(multialignment, output_dir, processing_time)

    if blocks_option:
        _create_blocks_file(multialignment, output_dir)

    for p in multialignment.poagraphs:
        poagraph_output_dir = _get_poagraph_output_dir(p, output_dir)
        _create_poagraph_sources_files(p, poagraph_output_dir)

        if consensus_option:
            _create_poagraph_consensus_files(p, poagraph_output_dir)

        if draw_poagraph_option:
            _create_poagraph_graph_files(p, poagraph_output_dir)


def _get_poagraph_output_dir(p, output_dir):
    return t.create_child_dir(output_dir, p.name)


def _create_common_files(multialignment, output_dir, processing_time):
    # wspolne pliki i informacje o multialignmencie (tu będą też bloki? grafy blokowe?)
    def copy_assets_dir():
        common_files_destination = t.join_path(output_dir, 'assets')
        common_files_source = t.join_path(WEBPAGE_DIR, 'assets')
        t.copy_dir(common_files_source, common_files_destination)

    def create_index_file():
        index_template_path = t.join_path(WEBPAGE_DIR, 'index.html')
        with open(index_template_path) as template_file:
            index_content = template_file.read()
            index_content = index_content.replace('info.js', str(multialignment.name) + """_info.js""")
            index_path = t.join_path(output_dir, 'index.html')

        with open(index_path, 'w') as output:
            output.write(index_content)

    def create_info_json():
        info_filename = t.join_path(output_dir, "info.json")
        info = {'name': multialignment.name,
                'running_time': processing_time,
                'sources_count': -1,
                'nodes_count': -1,
                'sequences_per_node': -1,
                'levels': [-1],
                'poagraphs': [src.name for src in multialignment.poagraphs]
                }
        with open(info_filename, 'w') as out_file:
            json.dump(info, fp=out_file, indent=4)

    copy_assets_dir()
    create_index_file()
    create_info_json()
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


def _create_blocks_file(multialignment, output_dir):
    nodes = []
    for node in multialignment.blocks_graph.nodes:
        nodes.append({'id': node})
        print(node)

    edges = []
    for (u, v) in multialignment.blocks_graph.edges:
        edges.append({"source": u,
                      "target": v,
                      "weight": multialignment.blocks_graph.edges[u, v]['weight'],
                      "active": multialignment.blocks_graph.edges[u, v]['active']})

    blocks = {"nodes": nodes, "edges": edges}

    blocks_filename = t.join_path(output_dir, "blocks.json")
    with open(blocks_filename, 'w') as out_file:
        json.dump(blocks, fp=out_file, indent=4)


def _create_poagraph_sources_files(poagraph, output_dir):
    sources_filename = t.join_path(output_dir, "sources.json")
    SourceEncoder.poagraph = poagraph
    with open(sources_filename, 'w') as out_file:
        json.dump(poagraph.sources, fp=out_file, cls=SourceEncoder, indent=4)


def _create_poagraph_consensus_files(poagraph, output_dir):
    consensuses_filename = t.join_path(output_dir, "consensuses.json")
    POAGraphRefEncoder.poagraph = poagraph
    with open(consensuses_filename, 'w') as out:
        json.dump(poagraph._poagraphrefs, fp=out, cls=POAGraphRefEncoder, indent=4)


def _create_poagraph_graph_files(poagraph, output_dir):
    # jeśli draw_poagraph_option - informacje o poagrafie
    def generate_blocks_graph(self, maf_file_path):
        # todo wykorzystać generowanie jsona
        # todo mergowanie bloków, które przechodzą w siebie po całości
        def get_nodes_list(blocks):
            #     def pies(id):
            #         return """{{
            #               data: {{
            #                 id: {0}
            #               }}
            #             }}""".format(id)
            #
            #     nodes_map = {}
            #     nodes_groups = [[]]
            #     for block in blocks:
            #         nextBlockids = set(block.srcID_to_next_blockID.values()) - set([None])
            #         if len(nextBlockids) == 1:
            #             iloop= 0
            #             for i, group in enumerate(nodes_groups):
            #                 iloop += 1
            #                 if block.ID in group:
            #                     nodes_groups[i].append(list(nextBlockids)[0])
            #             if iloop == len(nodes_groups):
            #                 nodes_groups.append([block.ID, list(nextBlockids)[0]])
            #         else:
            #             for r in nextBlockids:
            #                 nodes_groups.append([r])
            #     nodes_str_list = [pies(e) for e, group in enumerate(nodes_groups)]
            #     new_nodes_groups = []
            #     for g in nodes_groups:
            #         if g == []:
            #             continue
            #         else:
            #             new_nodes_groups.append(g)
            #     return new_nodes_groups, ",".join(nodes_str_list)
            nodes_list = []
            for block in blocks:
                nodes_list.append("""{{

                                                      data: {{
                                                        id: {0}
                                                      }}
                                                    }}""".format(block.ID))
                # nodes_list.append("""{{
                #                       data: {{
                #                         id: {0}
                #                       }},
                #                       position: {{ x: {1}, y: 0 }},
                #                     }}""".format(block.ID, block.ID *60 + 60))
            return [], ",".join(nodes_list)

        def get_edges_list(blocks, nodes_groups):
            def get_node_id(old_node_id):
                return old_node_id
                # for i, group in enumerate(nodes_groups):
                #     if old_node_id in group:
                #         return i

            def pies(i, src_dest, srcID):
                return """{{
                            data: {{
                                id: {0},
                                source: {1},
                                target: {2},
                                srcID: {3}
                                }},
                                classess: 'edge'
                                }}""".format(i, get_node_id(src_dest[0]), get_node_id(src_dest[1]), srcID)

            edges_list = {}
            edge_ID = len(blocks)
            for block in blocks:
                for srcID, next_blockID in block.srcID_to_next_blockID.items():
                    if next_blockID is None:
                        continue
                    if (block.ID, next_blockID) in edges_list.keys():
                        edges_list[(block.ID, next_blockID)].append(srcID)
                    else:
                        edge_ID += 1
                        edges_list[(block.ID, next_blockID)] = [srcID]

            edges_str_list = [pies(i + len(blocks), src_dest_srcID[0], src_dest_srcID[1]) for i, src_dest_srcID in
                              enumerate(edges_list.items())]

            return ",".join(edges_str_list)

        def get_blocks_data_as_json(blocks):
            nodes_groups, nodes_list = get_nodes_list(blocks)
            edges_list = get_edges_list(blocks, nodes_groups)
            return """var blocks = {{
                                        nodes: [{0}],
                                        edges: [{1}]}};""".format(nodes_list, edges_list)

        print("generating blocks graph")
        # prepare blocks
        self.name = self._get_multialignment_name(maf_file_path)
        self.output_dir = self._get_output_dir(maf_file_path)
        blocks = maf_reader.get_blocks(maf_file_path, self.name, self.output_dir)

        # prepare webpage template
        webpage_dir = t.get_real_path('../visualization_template')

        common_files_destination = t.join_path(self.output_dir, 'assets')
        common_files_source = t.join_path(webpage_dir, 'assets')
        t.copy_dir(common_files_source, common_files_destination)

        index_template_path = t.join_path(webpage_dir, 'index.html')
        with open(index_template_path) as template_file:
            index_content = template_file.read()
            index_content = index_content.replace('blocks.js', str(self.name) + """_blocks_data.js""")
            index_content = index_content.replace('info.js', str(self.name) + """_info.js""")

        index_path = t.join_path(self.output_dir, 'index.html')
        with open(index_path, 'w') as output:
            output.write(index_content)

        # prepare blocks data as json
        blocks_data_path = t.join_path(self.output_dir, str(self.name) + '_blocks_data.js')
        with open(blocks_data_path, 'w') as blocks_data_output:
            blocks_data_output.write(get_blocks_data_as_json(blocks))

    def convert_nodes_to_edges(poagraph):
        edges = []
        edge_id = -1
        for node in poagraph.nodes:
            for in_node in node.in_nodes:
                e = Edge(int(edge_id), int(in_node), node.ID, 1, -1, -1, "edge")
                edges.append(e)
                edge_id -= 1
            if node.aligned_to:
                e = Edge(int(edge_id), int(node.ID), int(node.aligned_to), 1, -1, -1, "edge aligned")
                edges.append(e)
                edge_id -=1
        for cons in poagraph.consensuses:
            cons_nodes = np.where(poagraph.nc[cons.ID, :]==True)[0]
            for i, n in enumerate(cons_nodes):
                if i == len(cons_nodes)-1:
                    break
                if n in poagraph.nodes[cons_nodes[i+1]].in_nodes:
                    e = Edge(id=int(edge_id),
                             source=int(n),
                             target=int(cons_nodes[i+1]),
                             weight=1,
                             consensus=cons.ID,
                             level=-1,
                             classes="edge consensus")
                    edges.append(e)
                    edge_id -= 1
            # sprawdzic przez jakie sources przechodzi
            # dla każdych nodow w tym sources zrobić krawędzie i odpowiednio je oznaczyć
        return edges

    # nodes
    nodes_filename = t.join_path(output_dir, "nodes.json")
    POAGraphNodeEncoder.poagraph = poagraph
    with open(nodes_filename, 'w') as out:
        json.dump(poagraph.nodes, fp=out, cls=POAGraphNodeEncoder, indent=4)

    # edges
    edges = convert_nodes_to_edges(poagraph)
    edges_filename = t.join_path(output_dir, "edges.json")
    POAGraphEdgeEncoder.poagraph = poagraph
    with open(edges_filename, 'w') as out:
        json.dump(edges, fp=out, cls=POAGraphEdgeEncoder, indent=4)