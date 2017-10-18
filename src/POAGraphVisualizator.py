import shutil
import os

import toolkit as t
from data_types import ebola as eb

class POAGraphVisualizator(object):
    def __init__(self, poagraph, output_dir, data_type='ebola'):
        self.output_dir = output_dir
        self.webpage_dir = t.get_real_path('../visualization_template')
        self.poagraph = poagraph
        self.data_type = data_type

    def generate(self, consensuses_comparison=False, graph_visualization=False):
        if consensuses_comparison:
            print("\tGenerating consensus comparison.")
            self.generate_consensus_comparison()
        elif graph_visualization:
            print("\tGenerating graph visualization.")
            self.generate_graph_visualization()
        else:
            return

        self.generate_common_files(consensuses_comparison, graph_visualization)


    def generate_consensus_comparison(self):
        sources_data_path = t.join_path(self.output_dir, str(self.poagraph.name)+'_sources_data.js')
        with open(sources_data_path, 'w') as sources_data_output:
            sources_data_output.write(self._get_sources_data_as_json())


    def generate_graph_visualization(self):
        pass


    def generate_common_files(self,consensuses_comparison, graph_visualization):
        common_files_destination = t.join_path(self.output_dir, 'assets')
        common_files_source = t.join_path(self.webpage_dir, 'assets')
        t.copy_dir(common_files_source, common_files_destination)

        index_template_path = t.join_path(self.webpage_dir, 'index.html')
        with open(index_template_path) as template_file:
            index_content = template_file.read()
            if consensuses_comparison:
                index_content = index_content.replace('sources_data.js', str(self.poagraph.name) + """_sources_data.js""")
            elif graph_visualization:
                index_content = index_content.replace('poagraph_data.js', str(self.poagraph.name) + """_poagraph_data.js""")


        index_path = t.join_path(self.output_dir, 'index.html')
        with open(index_path, 'w') as output:
            output.write(index_content)

    def _get_sources_data_as_json(self):
        sources_json = []
        consensuses_json = []
        poagraph_info = self._get_poagraph_info_json(self.poagraph)

        for source in self.poagraph.sources:
            sources_json.append(self._get_source_json(source))
        for consensus in self.poagraph.consensuses:
            consensuses_json.append(self._get_consensus_json(consensus))

        consensuses = "\n,".join(consensuses_json)
        sources = "\n,".join(sources_json)
        return """  var sources = {{ data: [\n{0}]}};\n
                    var consensuses = {{data: [\n{1}]}};\n
                    var poagraph = {2};\n""".format(sources, consensuses, poagraph_info)


    def _get_poagraph_info_json(self, poagraph):
        return """{{ nodes_count: {0},\tsources_count:{1},\tmean_sources_per_node: {2},\tconsensuses_count: {3}}}""".format \
                    (len(self.poagraph.nodes),
                    len(self.poagraph.sources),
                    self._get_mean_sources_per_node(),
                    len(self.poagraph.consensuses))


    def _get_mean_sources_per_node(self):
        def mean(numbers):
            return float(sum(numbers)) / max(len(numbers), 1)
        mean_sources_per_node = mean([node.sources_count for node in self.poagraph.nodes])
        return round(mean_sources_per_node, 2)

    def _get_source_json(self, source):
        if self.data_type is 'ebola':
            source_name = eb.extract_ebola_code_part(source.name, 0)
            source_alt_name = eb.extract_ebola_code_part(source.title, 1)
            source_group_name = eb.get_ebola_group_name(source_alt_name)
        else:
            source_name = source.name
            source_alt_name = '-'
            source_group_name = '-'

        return """{{ id: {0},\tname:'{1}',\ttitle: '{2}',\t group_name: '{3}',\tbundle_ID: {4},\tweight: {5},\tlength: {6}}}""".format(
                                                    source.ID,
                                                    source_name,
                                                    source_alt_name,
                                                    source_group_name,
                                                    source.consensusID,
                                                    source.weight,
                                                    len(source.nodes_IDs))

    def _get_consensus_json(self, consensus):
        return """{{ id: {0},\tname:'{1}',\ttitle: '{2}',\tsources_compatibility: {3},\tlength: {4} }}""".format\
                                                (consensus.ID,
                                                 consensus.name,
                                                 consensus.title,
                                                 consensus.compatibility_to_sources,
                                                 len(consensus.nodes_IDs))

    # def _get_poagraph_data_as_json(self, poagraph):
    #     poagraph_nodes_data = []
    #     poagraph_edges_data = []
    #     x_pos = 0
    #     edge_id = -1
    #     processed_nodes = []
    #     processed_edges = []
    #     srcs_count = len(poagraph.sources)
    #     consensus_ids = poagraph.consensuses.keys()
    #     processed_consensuses_edges = {consensus_ID: [] for consensus_ID in consensus_ids}
    #     consensus_indices = {consensus.ID : consensus.ID + len(poagraph.sources) for consensus in poagraph.consensuses.values()}
    #     for node_id, node in poagraph.nodedict.items():
    #         if node_id in processed_nodes:
    #             processed_nodes.remove(node_id)
    #             continue
    #         new_nodes_data, x_pos = self._get_node_data_as_json(node, poagraph, processed_nodes, x_pos, srcs_count, consensus_indices)
    #         poagraph_nodes_data += new_nodes_data
    #
    #     for node_id, node in poagraph.nodedict.items():
    #         new_edges_data, edge_id = self._get_edges_as_json(node, poagraph, processed_edges, processed_consensuses_edges, edge_id, srcs_count, consensus_ids)
    #         poagraph_edges_data += new_edges_data
    #
    #     for i, consensus_ID in enumerate(consensus_ids):
    #         new_edges_data, edge_id = self._get_consensus_edges_as_json(consensus_ID, i, poagraph, edge_id, srcs_count, consensus_indices)
    #         poagraph_edges_data += new_edges_data
    #
    #
    #     nodes = ',\n'.join(poagraph_nodes_data)
    #     edges = ',\n'.join(poagraph_edges_data)
    #     return """var PoaGraph_elements = {{
    #                         nodes: [\n{0}],
    #                         edges: [\n{1}]
    #                         }};""".format(nodes, edges)
    #
    #
    #
    # # nodes json generation
    #
    # def _get_node_data_as_json(self, node, poagraph, processed_nodes, x_pos, srcs_count, consensus_indices):
    #     json = []
    #     x_pos, y_pos = self._calc_node_pos(node, poagraph, x_pos)
    #     if node.alignedTo:
    #         aligned_count = len(node.alignedTo)
    #         if aligned_count == 3:
    #             y_positions = [-30, -15, 15, 30]
    #         elif aligned_count == 2:
    #             y_positions = [-15, 0, 15]
    #         elif aligned_count == 1:
    #             y_positions = [-15, 15]
    #         json.append(
    #             self._get_node_json(node, self._get_node_weight(node.ID, poagraph, srcs_count), x_pos, y_positions[0],consensus_indices))
    #         #processed_nodes.append(node.ID)
    #         for i, aligned_id in enumerate(node.alignedTo):
    #             aligned_node = poagraph.nodedict[aligned_id]
    #             json.append(
    #                 self._get_node_json(aligned_node, self._get_node_weight(aligned_node.ID, poagraph, srcs_count),
    #                                     x_pos, y_positions[i + 1], consensus_indices))
    #             processed_nodes.append(aligned_node.ID)
    #     json.append(self._get_node_json(node, self._get_node_weight(node.ID, poagraph, srcs_count), x_pos, y_pos,consensus_indices))
    #
    #     return (json, x_pos)
    #
    #
    # def _calc_node_pos(self, node, poagraph, current_x_pos):
    #     if len(node.sourceSequenceIds) == len(poagraph.sources):
    #         return (current_x_pos+35, 0)
    #     y = (node.sourceSequenceIds[0] - (-10)) / (10 - (-10))
    #     return (current_x_pos + 60, y)
    #
    #
    # def _get_node_weight(self, node_ID, poagraph, srcs_count_no_hb):
    #     node_src_count_no_hb = 0
    #     for src_ID in poagraph.nodedict[node_ID].sourceSequenceIds:
    #         node_src_count_no_hb += 1
    #         # if not poagraph.sources[src_ID].is_consensus:
    #         #     node_src_count_no_hb += 1
    #     return node_src_count_no_hb/srcs_count_no_hb
    #
    #
    # def _get_node_json(self, node, node_weight, x_pos, y_pos, consensus_indices):
    #     json = """{{data: {{ id: {0},    nucleobase: '{1}', source: {2}, weight: {3} }},
    #                 position: {{ x: {4}, y: {5} }}
    #                 }}""".format(node.ID, node.base, self._get_sources_as_table(node, consensus_indices), node_weight, x_pos, y_pos)
    #     return json
    #
    #
    # def _get_sources_as_table(self, node, consensus_indices):
    #     sources = '['
    #     for source in node.sourceSequenceIds:
    #         sources += (str(source) + ',')
    #     #TODO dwie poniższe linie dodałam
    #     for consensusID in node.consensusIds:
    #         sources += (str(consensus_indices[consensusID])+',')
    #     sources = sources[:-1] + ']'
    #     return sources
    #
    #
    # # edges json generation
    #
    # def _get_edges_as_json(self, node, poagraph, processed_edges, processed_consensus_edges, edge_id, srcs_count, consensuses_ids):
    #     json = []
    #     edge_type = 0
    #     # normal edges entering this node
    #     for in_node_ID in node.inNodes:
    #         json.append(self._get_edge_json(src=in_node_ID, target=node.ID, consensus_id=-1, edge_type=edge_type, edge_id=edge_id, edge_weight=self._get_edge_weight(node.ID, in_node_ID, poagraph, srcs_count)))
    #         edge_id -= 1
    #
    #     # edges to aligned nodes
    #     if node.alignedTo:
    #         nearest_aligned_node = min(node.alignedTo)
    #         # for aligned_node_ID in node.alignedTo[0:1]:
    #         if (node.ID, nearest_aligned_node) in processed_edges:
    #             processed_edges.remove((node.ID, nearest_aligned_node))
    #         elif (nearest_aligned_node, node.ID) in processed_edges:
    #             processed_edges.remove((nearest_aligned_node, node.ID))
    #         else:
    #             json.append(self._get_edge_json(src=nearest_aligned_node, target=node.ID, consensus_id=-1, edge_type=1, edge_id=edge_id, edge_weight=1))
    #             edge_id -= 1
    #             processed_edges.append((node.ID, nearest_aligned_node))
    #
    #     return (json, edge_id)
    #
    #
    # def _get_edge_weight(self, src_node_ID, target_node_ID, poagraph, srcs_count):
    #     source_node = poagraph.nodedict[src_node_ID]
    #     target_node = poagraph.nodedict[target_node_ID]
    #     common_source_sequences_count = 0
    #     for source_sequence_id in source_node.sourceSequenceIds:
    #         #if not poagraph.sources[source_sequence_id].is_consensus and source_sequence_id in target_node.sourceSequenceIds:
    #         if source_sequence_id in target_node.sourceSequenceIds:
    #             common_source_sequences_count += 1
    #
    #     return common_source_sequences_count/srcs_count
    #
    # def _get_edge_json(self, src, target, consensus_id, edge_type, edge_id, edge_weight):
    #     #0 type - normalna krawedz,
    #     #1 type - krawedz aligned - wtedy klasa aligned
    #     #2 type - consensus
    #     css_class = """'edge'"""
    #     if edge_type == 1:
    #         css_class = """'edge aligned'"""
    #     elif edge_type == 2:
    #         css_class = """'edge consensus'"""
    #     json = """{{data: {{ id: {0},    source: {1}, target: {2}, consensus: {3},    weight: {4} }},
    #                 classes: {5}
    #                 }}""".format(edge_id,
    #                              src,
    #                              target,
    #                              consensus_id,
    #                              edge_weight,
    #                              css_class)
    #     return json
    #
    #



    #

    #
    #

    #
    #

    #
    #
    # def _get_index_html(self, poagraph_id):
    #     index_template_path = ''.join([self.webpage_dir, '/index.html'])
    #     with open(index_template_path) as template_file:
    #         index_content = template_file.read()
    #         index_content = index_content.replace('poagraph_data.js', str(poagraph_id) + """_poagraph_data.js""")
    #         index_content = index_content.replace('sources_data.js', str(poagraph_id) + """_sources_data.js""")
    #     return index_content
    #
    #
    # def _get_consensus_edges_as_json(self, consensus_big_ID, consensus_small_ID, poagraph, edge_id, srcs_count, consensus_indices):
    #     #consensus_big_ID = consensus_indices[consensus_small_ID]
    #     json = []
    #     consensus_nodes_IDs = []
    #     for node_ID, node in poagraph.nodedict.items():
    #         if consensus_small_ID in node.consensusIds:
    #             consensus_nodes_IDs.append(node_ID)
    #         # if consensus_big_ID in node.sourceSequenceIds:
    #         #     consensus_nodes_IDs.append(node_ID)
    #     for i, node_ID in enumerate(consensus_nodes_IDs):
    #         if i < len(consensus_nodes_IDs)-1:
    #             json.append(self._get_edge_json(src=node_ID, target=consensus_nodes_IDs[i+1], consensus_id=consensus_small_ID, edge_type=2,
    #                                         edge_id=edge_id,
    #                                         edge_weight=self._get_edge_weight(node_ID, consensus_nodes_IDs[i+1], poagraph, srcs_count)))
    #             edge_id -= 1
    #     return (json, edge_id)
    #
    #
    #
