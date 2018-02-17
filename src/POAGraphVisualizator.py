# import toolkit as t
# from data_types import ebola as eb
#
# class POAGraphVisualizator(object):
#     def __init__(self, poagraph, output_dir, data_type='ebola'):
#         self.output_dir = output_dir
#         self.webpage_dir = t.get_real_path('../visualization_template')
#         self.poagraph = poagraph
#         self.data_type = data_type
#
#     def generate(self, consensuses_comparison=False, graph_visualization=False, processing_time = '', tresholds='', consensus_algorithm=False):
#         if consensuses_comparison:
#             print("\tGenerating consensuses info.")
#             self.generate_consensus_comparison()
#         if graph_visualization:
#             print("\tGenerating graph visualization.")
#             self.generate_graph_visualization()
#
#         print("\tGenerating info file.")
#         self.generate_info_file(processing_time, tresholds, consensus_algorithm)
#
#         self.generate_common_files(consensuses_comparison, graph_visualization)
#
#     def generate_consensus_comparison(self):
#         for c in self.poagraph.consensuses:
#             if c.parent_consensus != -1:
#                 self.poagraph.consensuses[c.parent_consensus].children.append(c)
#
#         sources_data_path = t.join_path(self.output_dir, str(self.poagraph.name)+'_sources_data.js')
#         with open(sources_data_path, 'w') as sources_data_output:
#             sources_data_output.write(self._get_sources_data_as_json())
#
#         consensuses_data_path = t.join_path(self.output_dir, str(self.poagraph.name)+'_consensuses_data.js')
#         with open(consensuses_data_path, 'w') as consensuses_data_output:
#             consensuses_data_output.write(self._get_consensuses_data_as_json())
#
#     def generate_graph_visualization(self):
#         poagraph_data_path = t.join_path(self.output_dir, str(self.poagraph.name) + '_poagraph_data.js')
#         with open(poagraph_data_path, 'w') as poagraph_data_output:
#             poagraph_data_output.write(self._get_poagraph_data_as_json())
#
#     def generate_info_file(self, processing_time, tresholds, consensus_algorithm):
#         info_path = t.join_path(self.output_dir, str(self.poagraph.name) + '_info.js')
#         with open(info_path, 'w') as info_output:
#             info_output.write(self._get_info_as_json(processing_time, tresholds, consensus_algorithm))
#
    # def generate_common_files(self, consensuses_comparison, graph_visualization):
    #     common_files_destination = t.join_path(self.output_dir, 'assets')
    #     common_files_source = t.join_path(self.webpage_dir, 'assets')
    #     t.copy_dir(common_files_source, common_files_destination)
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
#
#     def _get_sources_data_as_json(self):
#         sources_json = []
#         for source in self.poagraph.sources:
#             sources_json.append(self._get_source_json(source))
#
#         sources = "\n,".join(sources_json)
#         return """  var sources = {{ data: [\n{0}]}};""".format(sources)
#
#     def _get_consensuses_data_as_json(self):
#
#         consensues_json = []
#         for consensus in self.poagraph.consensuses:
#             consensues_json.append(self._get_consensus_json(consensus))
#         consensuses_json_table = ','.join(consensues_json)
#
#         return """  var consensuses = {{ c : [ {0}]}};""".format(consensuses_json_table)
#         # consensuses_levels = set([consensus.level for consensus in self.poagraph.consensuses])
#
#         # consensuses_levels_json = []
#         # for level in consensuses_levels:
#         #     current_level_consensuses = [consensus for consensus in self.poagraph.consensuses if consensus.level == level]
#         #     consensuses_levels_json.append(self._get_consensuses_level_json(current_level_consensuses))
#
#         # consensuses = "\n,".join(consensuses_levels_json)
#         # return """  var consensuses = {{ c : [ {0}]}};""".format(consensuses)
#
#     def _get_consensuses_level_json(self, current_level_consensuses):
#         def get_consensuses_json(current_level_consensuses):
#             consensues_json = []
#             for consensus in current_level_consensuses:
#                 consensues_json.append(self._get_consensus_json(consensus))
#             return ','.join(consensues_json)
#
#         return """{{ treshold: {0},\tdata: [{1}]}}""".format(str(current_level_consensuses[0].level),
#                                                          get_consensuses_json(current_level_consensuses))
#
#     def _get_info_as_json(self, processing_time, tresholds, consensus_algorithm):
#         return """var info =
#                 {{ name: '{0}', \nconsensus_algorithm: {1}, \nrunning_time:'{2}', \nsources_count: {3}, \nnodes_count: {4}, \nsequences_per_node: {5}, \nlevels: {6}}}; \nvar treeData = {7};""".format \
#                     (self.poagraph.name,
#                      consensus_algorithm if consensus_algorithm else '',
#                      processing_time,
#                      len(self.poagraph.sources),
#                      len(self.poagraph.nodes),
#                      self._get_mean_sources_per_node(),
#                      tresholds,
#                      self._get_consensuses_as_tree(tresholds)
#                      )
#
#     def _get_mean_sources_per_node(self):
#         def mean(numbers):
#             return float(sum(numbers)) / max(len(numbers), 1)
#
#         mean_sources_per_node = mean([len(node.sources) for node in self.poagraph.nodes])
#         return round(mean_sources_per_node, 2)
#
#     def _get_source_json(self, source):
#         if self.data_type is 'ebola':
#             source_name = eb.extract_ebola_code_part(source.name, 0)
#             source_alt_name = eb.get_official_ebola_name(eb.extract_ebola_code_part(source.title, 1))
#             source_group_name = eb.get_ebola_group_name(eb.extract_ebola_code_part(source.title, 1))
#         else:
#             source_name = source.name
#             source_alt_name = '-'
#             source_group_name = '-'
#
#         return """{{ id: {0},\tname:'{1}',\ttitle: '{2}',\t group: '{3}',\tbundle_ID: {4},\tweight: {5},\tfirst_node_id: {6},\tlength: {7}}}""".format(
#                                                     source.currentID,
#                                                     source_name,
#                                                     source_alt_name,
#                                                     source_group_name,
#                                                     str(source.consensuses),
#                                                     source.weight,
#                                                     source.nodes_IDs[0],
#                                                     len(source.nodes_IDs))
#
#     def _get_consensus_json(self, consensus):
#         consensus_children = [c.currentID for c in consensus.children]
#         return """{{ id: {0},\tname:'{1}',\ttitle: '{2}',\tsources_compatibility: {3},\tlength: {4},\tsources: {5},\tparent: {6},\tlevel:{7}, \tchildren:{8}}}""".format\
#                                                 (consensus.currentID,
#                                                  consensus.name,
#                                                  consensus.title,
#                                                  consensus.compatibility_to_sources,
#                                                  len(consensus.nodes_IDs),
#                                                  consensus.sources_IDs,
#                                                  consensus.parent_consensus,
#                                                  consensus.level,
#                                                  consensus_children)
#
#     def _get_consensuses_as_tree(self, tresholds):
#         # def get_leaves(sources_ids, parent):
#         #     l =  """[""";
#         #     kk = []
#         #     for i in sources_ids:
#         #         kk.append("""{{val: 1, name: '{0}', \nparent: '{1}'}}""".format("_".join(["source", str(i)]), parent))
#         #     kk = ", ".join(kk)
#         #     l += kk
#         #     l+= """]"""
#         #     return l
#
#         def get_consensus_node(consensus):
#             if not consensus.children:
#                 return """{{name: '{0}', \nparent: '{1}', \nval: {2}, \nsources: {3}}}""".format(consensus.name,
#                                                                                                   consensus.parent_consensus,
#                                                                                                   consensus.level,
#                                                                                                   #get_leaves(consensus.sources_IDs, consensus.currentID)),
#                                                                                                   consensus.sources_IDs)
#
#
#             children = []
#             for child in consensus.children:
#                 children.append(get_consensus_node(child))
#
#             children = ",".join(children)
#             return """{{name: '{0}', \nparent:'{1}', \nval: {2}, \nchildren:[{3}]}}""".format(consensus.name, str(consensus.parent_consensus), consensus.level, children)
#
#         # for c in self.poagraph.consensuses:
#         #     if c.parent_consensus != -1:
#         #         self.poagraph.consensuses[c.parent_consensus].children.append(c)
#
#         ccc = []
#         for c in self.poagraph.consensuses:
#             if c.parent_consensus == -1:
#                 ccc.append(get_consensus_node(c))
#
#         ccc_json = ",".join(ccc)
#
#         return """[{{name: 'All sequences', \nparent: 'null', \nval: 0, \nchildren: [{0}]}}]""".format(ccc_json)
#         # return """[
#         #           {
#         #             "name": "Top Level",
#         #             "parent": "null",
#         #             "children": [
#         #               {
#         #                 "name": "Level 2: A",
#         #                 "parent": "Top Level",
#         #                 "children": [
#         #                   {
#         #                     "name": "Son of A",
#         #                     "parent": "Level 2: A"
#         #                   },
#         #                   {
#         #                     "name": "Daughter of A",
#         #                     "parent": "Level 2: A"
#         #                   }
#         #                 ]
#         #               },
#         #               {
#         #                 "name": "Level 2: B",
#         #                 "parent": "Top Level"
#         #               }
#         #             ]
#         #           }
#         #         ]"""
#
#     def _get_poagraph_data_as_json(self):
#         poagraph_nodes_data = []
#         poagraph_edges_data = []
#         x_pos = 0
#         edge_id = -1
#         processed_nodes = []
#         processed_edges = []
#         srcs_count = len(self.poagraph.sources)
#         consensus_ids = [consensus.currentID for consensus in self.poagraph.consensuses if consensus.level == 0]
#         processed_consensuses_edges = {consensus_ID: [] for consensus_ID in consensus_ids}
#         consensus_indices = {consensus.currentID : consensus.currentID + len(self.poagraph.sources) for consensus in self.poagraph.consensuses if consensus.level == 0}
#         for node in self.poagraph.nodes:
#             if node.currentID in processed_nodes:
#                 processed_nodes.remove(node.currentID)
#                 continue
#             new_nodes_data, x_pos = self._get_node_data_as_json(node, processed_nodes, x_pos, srcs_count, consensus_indices)
#             poagraph_nodes_data += new_nodes_data
#
#         for node in self.poagraph.nodes:
#             new_edges_data, edge_id = self._get_edges_as_json(node, processed_edges, processed_consensuses_edges, edge_id, srcs_count, consensus_ids)
#             poagraph_edges_data += new_edges_data
#
#         for i, consensus_ID in enumerate(consensus_ids):
#             new_edges_data, edge_id = self._get_consensus_edges_as_json(consensus_ID, i, self.poagraph, edge_id, srcs_count, consensus_indices)
#             poagraph_edges_data += new_edges_data
#
#
#         nodes = ',\n'.join(poagraph_nodes_data)
#         edges = ',\n'.join(poagraph_edges_data)
#         return """var poagraph = {{
#                             nodes: [\n{0}],
#                             edges: [\n{1}]
#                             }};""".format(nodes, edges)
#
#     def _get_node_data_as_json(self, node, processed_nodes, x_pos, srcs_count, consensus_indices):
#         json = []
#         x_pos, y_pos = self._calc_node_pos(node, self.poagraph, x_pos)
#         if node.aligned_to:
#             aligned_count = len(node.aligned_to)
#             if aligned_count == 3:
#                 y_positions = [-30, -15, 15, 30]
#             elif aligned_count == 2:
#                 y_positions = [-15, 0, 15]
#             elif aligned_count == 1:
#                 y_positions = [-15, 15]
#             json.append(
#                 self._get_node_json(node, self._get_node_weight(node.currentID, self.poagraph, srcs_count), x_pos, y_positions[0], consensus_indices))
#             #processed_nodes.append(node.ID)
#             for i, aligned_id in enumerate(node.aligned_to):
#                 aligned_node = self.poagraph.nodes[aligned_id]
#                 json.append(
#                     self._get_node_json(aligned_node, self._get_node_weight(aligned_node.currentID, self.poagraph, srcs_count),
#                                         x_pos, y_positions[i + 1], consensus_indices))
#                 processed_nodes.append(aligned_node.currentID)
#         json.append(self._get_node_json(node, self._get_node_weight(node.currentID, self.poagraph, srcs_count), x_pos, y_pos, consensus_indices))
#
#         return (json, x_pos)
#
#     def _calc_node_pos(self, node, poagraph, current_x_pos):
#         if len(node.sources) == len(poagraph.sources):
#             return (current_x_pos+35, 0)
#         node_sources = sorted(node.sources)
#         y = (node_sources[0] - (-10)) / (10 - (-10))
#         return (current_x_pos + 60, y)
#
#     def _get_node_weight(self, node_ID, poagraph, srcs_count_no_hb):
#         node_src_count_no_hb = 0
#         for src_ID in poagraph.nodes[node_ID].sources:
#             node_src_count_no_hb += 1
#             # if not poagraph.sources[src_ID].is_consensus:
#             #     node_src_count_no_hb += 1
#         return node_src_count_no_hb/srcs_count_no_hb
#
#     def _get_node_json(self, node, node_weight, x_pos, y_pos, consensus_indices):
#         json = """{{data: {{ id: {0},    nucleobase: '{1}', source: {2}, weight: {3} }},
#                     position: {{ x: {4}, y: {5} }}
#                     }}""".format(node.currentID, node.base, self._get_sources_as_table(node, consensus_indices), node_weight, x_pos, y_pos)
#         return json
#
#     def _get_sources_as_table(self, node, consensus_indices):
#         sources = '['
#         node_sources = sorted(node.sources)
#         for source in node_sources:
#             sources += (str(source) + ',')
#         #TODO dwie poniższe linie dodałam
#         current_node_consensuses = [consensus.currentID for consensus in self.poagraph.consensuses if node.currentID in consensus.nodes_IDs and consensus.level == 0]
#         for consensusID in current_node_consensuses:
#             sources += (str(consensus_indices[consensusID])+',')
#         sources = sources[:-1] + ']'
#         return sources
#
#     # edges json generation
#
#     def _get_edges_as_json(self, node, processed_edges, processed_consensus_edges, edge_id, srcs_count, consensuses_ids):
#         json = []
#         edge_type = 0
#         # normal edges entering this node
#         for in_node_ID in node.in_nodes:
#             json.append(self._get_edge_json(src=in_node_ID, target=node.currentID, consensus_id=-1, edge_type=edge_type, edge_id=edge_id, edge_weight=self._get_edge_weight(node.currentID, in_node_ID, self.poagraph, srcs_count)))
#             edge_id -= 1
#
#         # edges to aligned nodes
#         if node.aligned_to:
#             nearest_aligned_node = min(node.aligned_to)
#             # for aligned_node_ID in node.alignedTo[0:1]:
#             if (node.currentID, nearest_aligned_node) in processed_edges:
#                 processed_edges.remove((node.currentID, nearest_aligned_node))
#             elif (nearest_aligned_node, node.currentID) in processed_edges:
#                 processed_edges.remove((nearest_aligned_node, node.currentID))
#             else:
#                 json.append(self._get_edge_json(src=nearest_aligned_node, target=node.currentID, consensus_id=-1, edge_type=1, edge_id=edge_id, edge_weight=1))
#                 edge_id -= 1
#                 processed_edges.append((node.currentID, nearest_aligned_node))
#
#         return (json, edge_id)
#
#     def _get_edge_weight(self, src_node_ID, target_node_ID, poagraph, srcs_count):
#         source_node = poagraph.nodes[src_node_ID]
#         target_node = poagraph.nodes[target_node_ID]
#         common_source_sequences_count = 0
#         source_node_sources = sorted(source_node.sources)
#         target_node_sources = sorted(target_node.sources)
#         for source_sequence_id in source_node_sources:
#             #if not poagraph.sources[source_sequence_id].is_consensus and source_sequence_id in target_node.sourceSequenceIds:
#
#             if source_sequence_id in target_node_sources:
#                 common_source_sequences_count += 1
#
#         return common_source_sequences_count/srcs_count
#
#     def _get_edge_json(self, src, target, consensus_id, edge_type, edge_id, edge_weight):
#         #0 type - normalna krawedz,
#         #1 type - krawedz aligned - wtedy klasa aligned
#         #2 type - consensus
#         css_class = """'edge'"""
#         if edge_type == 1:
#             css_class = """'edge aligned'"""
#         elif edge_type == 2:
#             css_class = """'edge consensus'"""
#         json = """{{data: {{ id: {0},    source: {1}, target: {2}, consensus: {3},    weight: {4} }},
#                     classes: {5}
#                     }}""".format(edge_id,
#                                  src,
#                                  target,
#                                  consensus_id,
#                                  edge_weight,
#                                  css_class)
#         return json
#
#     def _get_consensus_edges_as_json(self, consensus_big_ID, consensus_small_ID, poagraph, edge_id, srcs_count, consensus_indices):
#         #consensus_big_ID = consensus_indices[consensus_small_ID]
#         json = []
#         consensus_nodes_IDs = []
#         for node in poagraph.nodes:
#             node_consensuses = [consensus.currentID for consensus in self.poagraph.consensuses if node.currentID in consensus.nodes_IDs and consensus.level == 0]
#             if consensus_small_ID in node_consensuses:
#                 consensus_nodes_IDs.append(node.currentID)
#             # if consensus_big_ID in node.sourceSequenceIds:
#             #     consensus_nodes_IDs.append(node_ID)
#         for i, node_ID in enumerate(consensus_nodes_IDs):
#             if i < len(consensus_nodes_IDs)-1:
#                 json.append(self._get_edge_json(src=node_ID, target=consensus_nodes_IDs[i+1], consensus_id=consensus_small_ID, edge_type=2,
#                                             edge_id=edge_id,
#                                             edge_weight=self._get_edge_weight(node_ID, consensus_nodes_IDs[i+1], poagraph, srcs_count)))
#                 edge_id -= 1
#         return (json, edge_id)