import toolkit as t
from data_types import ebola as eb

class POAGraphVisualizator(object):
    def __init__(self, poagraph, output_dir, data_type='ebola'):
        self.output_dir = output_dir
        self.webpage_dir = t.get_real_path('../visualization_template_new')
        self.poagraph = poagraph
        self.data_type = data_type

    def generate(self, consensuses_comparison=False, graph_visualization=False, processing_time = '', tresholds='', consensus_algorithm=False):
        if consensuses_comparison:
            print("\tGenerating consensuses info.")
            self.generate_consensus_comparison()
        elif graph_visualization:
            print("\tGenerating graph visualization.")
            self.generate_graph_visualization()

        print("\tGenerating info file.")
        self.generate_info_file(processing_time, tresholds, consensus_algorithm)

        self.generate_common_files(consensuses_comparison, graph_visualization)


    def generate_consensus_comparison(self):
        sources_data_path = t.join_path(self.output_dir, str(self.poagraph.name)+'_sources_data.js')
        with open(sources_data_path, 'w') as sources_data_output:
            sources_data_output.write(self._get_sources_data_as_json())

        consensuses_data_path = t.join_path(self.output_dir, str(self.poagraph.name)+'_consensuses_data.js')
        with open(consensuses_data_path, 'w') as consensuses_data_output:
            consensuses_data_output.write(self._get_consensuses_data_as_json())

    def generate_graph_visualization(self):
        pass

    def generate_info_file(self, processing_time, tresholds, consensus_algorithm):
        info_path = t.join_path(self.output_dir, str(self.poagraph.name) + '_info.js')
        with open(info_path, 'w') as info_output:
            info_output.write(self._get_info_as_json(processing_time, tresholds, consensus_algorithm))

    def generate_common_files(self, consensuses_comparison, graph_visualization):
        common_files_destination = t.join_path(self.output_dir, 'assets')
        common_files_source = t.join_path(self.webpage_dir, 'assets')
        t.copy_dir(common_files_source, common_files_destination)

        index_template_path = t.join_path(self.webpage_dir, 'index.html')
        with open(index_template_path) as template_file:
            index_content = template_file.read()
            if consensuses_comparison:
                index_content = index_content.replace('sources.js', str(self.poagraph.name) + """_sources_data.js""")
                index_content = index_content.replace('consensuses.js', str(self.poagraph.name) + """_consensuses_data.js""")
            elif graph_visualization:
                index_content = index_content.replace('poagraph.js', str(self.poagraph.name) + """_poagraph_data.js""")
            index_content = index_content.replace('info.js', str(self.poagraph.name) + """_info.js""")

        index_path = t.join_path(self.output_dir, 'index.html')
        with open(index_path, 'w') as output:
            output.write(index_content)

    def _get_sources_data_as_json(self):
        sources_json = []
        for source in self.poagraph.sources:
            sources_json.append(self._get_source_json(source))

        sources = "\n,".join(sources_json)
        return """  var sources = {{ data: [\n{0}]}};""".format(sources)

    def _get_consensuses_data_as_json(self):
        consensuses_levels = set([consensus.level for consensus in self.poagraph.consensuses])

        consensuses_levels_json = []
        for level in consensuses_levels:
            current_level_consensuses = [consensus for consensus in self.poagraph.consensuses if consensus.level == level]
            consensuses_levels_json.append(self._get_consensuses_level_json(current_level_consensuses))

        consensuses = "\n,".join(consensuses_levels_json)
        return """  var consensuses = {{ c : [ {0}]}};""".format(consensuses)

    def _get_consensuses_level_json(self, current_level_consensuses):
        def get_consensuses_json(current_level_consensuses):
            consensues_json = []
            for consensus in current_level_consensuses:
                consensues_json.append(self._get_consensus_json(consensus))
            return ','.join(consensues_json)

        return """{{ treshold: {0},\tdata: [{1}]}}""".format(str(current_level_consensuses[0].level),
                                                         get_consensuses_json(current_level_consensuses))

    def _get_info_as_json(self, processing_time, tresholds, consensus_algorithm):
        return """var info = 
                {{ name: '{0}', \nconsensus_algorithm: {1}, \nrunning_time:'{2}', \nsources_count: {3}, \nnodes_count: {4}, \nsequences_per_node: {5}, \nlevels: {6}}}; \nvar newick = '{7}';""".format \
                    (self.poagraph.name,
                     consensus_algorithm if consensus_algorithm else '',
                     processing_time,
                     len(self.poagraph.sources),
                     len(self.poagraph.nodes),
                     self._get_mean_sources_per_node(),
                     tresholds,
                     self._get_consensuses_as_newick()
                     )

    def _get_mean_sources_per_node(self):
        def mean(numbers):
            return float(sum(numbers)) / max(len(numbers), 1)
        mean_sources_per_node = mean([len(node.sources) for node in self.poagraph.nodes])
        return round(mean_sources_per_node, 2)

    def _get_source_json(self, source):
        if self.data_type is 'ebola':
            source_name = eb.extract_ebola_code_part(source.name, 0)
            source_alt_name = eb.get_official_ebola_name(eb.extract_ebola_code_part(source.title, 1))
            source_group_name = eb.get_ebola_group_name(eb.extract_ebola_code_part(source.title, 1))
        else:
            source_name = source.name
            source_alt_name = '-'
            source_group_name = '-'

        return """{{ id: {0},\tname:'{1}',\ttitle: '{2}',\t group: '{3}',\tbundle_ID: {4},\tweight: {5},\tfirst_node_id: {6},\tlength: {7}}}""".format(
                                                    source.currentID,
                                                    source_name,
                                                    source_alt_name,
                                                    source_group_name,
                                                    str(source.consensuses),
                                                    source.weight,
                                                    source.nodes_IDs[0],
                                                    len(source.nodes_IDs))

    def _get_consensus_json(self, consensus):
        return """{{ id: {0},\tname:'{1}',\ttitle: '{2}',\tsources_compatibility: {3},\tlength: {4} }}""".format\
                                                (consensus.currentID,
                                                 consensus.name,
                                                 consensus.title,
                                                 consensus.compatibility_to_sources,
                                                 len(consensus.nodes_IDs))

    def _get_consensuses_as_newick(self):
        
        return "(((EELA:0.150276,CONGERA:0.213019):0.230956,(EELB:0.263487,CONGERB:0.202633):0.246917):0.094785,((CAVEFISH:0.451027,(GOLDFISH:0.340495,ZEBRAFISH:0.390163):0.220565):0.067778,((((((NSAM:0.008113,NARG:0.014065):0.052991,SPUN:0.061003,(SMIC:0.027806,SDIA:0.015298,SXAN:0.046873):0.046977):0.009822,(NAUR:0.081298,(SSPI:0.023876,STIE:0.013652):0.058179):0.091775):0.073346,(MVIO:0.012271,MBER:0.039798):0.178835):0.147992,((BFNKILLIFISH:0.317455,(ONIL:0.029217,XCAU:0.084388):0.201166):0.055908,THORNYHEAD:0.252481):0.061905):0.157214,LAMPFISH:0.717196,((SCABBARDA:0.189684,SCABBARDB:0.362015):0.282263,((VIPERFISH:0.318217,BLACKDRAGON:0.109912):0.123642,LOOSEJAW:0.397100):0.287152):0.140663):0.206729):0.222485,(COELACANTH:0.558103,((CLAWEDFROG:0.441842,SALAMANDER:0.299607):0.135307,((CHAMELEON:0.771665,((PIGEON:0.150909,CHICKEN:0.172733):0.082163,ZEBRAFINCH:0.099172):0.272338):0.014055,((BOVINE:0.167569,DOLPHIN:0.157450):0.104783,ELEPHANT:0.166557):0.367205):0.050892):0.114731):0.295021)"
