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

