# import os
from subprocess import run
import toolkit as t
import maf_reader as maf_reader
import po_reader as po_reader
from POAGraphVisualizator import POAGraphVisualizator
from Sequence import Consensus
# from fasta_generators import generate_source_as_fasta_from_poagraph, generate_consensus_as_fasta_from_poagraph

class Multialignment(object):
    def __init__(self, data_type):
        self.name = None
        self.output_dir = None
        self.poagraphs = None
        self.data_type=data_type


    def build_multialignment_from_maf(self, maf_file_name, merge_option):
        print("Buliding multialignment from " + maf_file_name + "...")
        def _get_mutlialignment_name(input_file_name):
            return t.get_file_name_without_extension(input_file_name)

        def _get_output_dir(input_file_name):
            return t.create_next_sibling_dir(input_file_name, "converted")

        self.name = _get_mutlialignment_name(maf_file_name)
        self.output_dir = _get_output_dir(maf_file_name)
        self.poagraphs = maf_reader.parse_to_poagraphs(file_path = maf_file_name,
                                                       merge_option = merge_option,
                                                       multialignment_name = self.name,
                                                       output_dir = self.output_dir)

    def build_multialignment_from_po(self, po_file_name):
        print("Buliding multialignment from " + po_file_name + "...")
        def _get_mutlialignment_name(input_file_name):
            return t.get_file_name_without_extension(input_file_name)

        def _get_output_dir(input_file_name):
            return t.create_next_sibling_dir(input_file_name, "converted")

        self.name = _get_mutlialignment_name(po_file_name)
        self.output_dir = _get_output_dir(po_file_name)
        self.poagraphs = [po_reader.parse_to_poagraph(file_path = po_file_name,
                                                        output_dir = self.output_dir)]


    def generate_consensus(self, consensus_iterative, hbmin, min_comp):
        for i, p in enumerate(self.poagraphs):
            consensus_output_dir = t.create_child_dir(p.path, "consensus")
            if consensus_iterative:
                print('Generate consensuses (hbmin='+ str(hbmin) +', min_comp=' + str(min_comp) + ') iteratively...')
                self._run_iterative_consensus_generation(consensus_output_dir, p, hbmin, min_comp)
            else:
                print('Generate consensuses (hbmin='+ str(hbmin) +') in one iteration...')
                new_poagraph = self._run_single_consensus_generation(consensus_output_dir, p, hbmin)
                new_poagraph.path = p.path
                self.poagraphs[i] = new_poagraph


    def _run_single_consensus_generation(self, consensus_output_dir, poagraph, hbmin, consensus_name = "consensus"):
        poagraph_as_po = poagraph.generate_po()

        file_name = t.join_path(consensus_output_dir, consensus_name)
        hb_file_name = t.change_file_extension(file_name,  '.hb')
        with open(file_name, 'w') as output_po_file:
            output_po_file.write(poagraph_as_po)
        run(['../bin/poa', '-read_msa', file_name, '-hb', '-po', hb_file_name, '../bin/blosum80.mat', '-hbmin',
             str(hbmin)])

        new_poagraph = po_reader.parse_to_poagraph(hb_file_name, poagraph.path)

        if all([src.consensusID == -1 for src in new_poagraph.sources]):
            raise NoConsensusFound()

        new_poagraph.calculate_compatibility_to_consensuses()
        return new_poagraph

    def _run_iterative_consensus_generation(self, consensus_output_dir, poagraph, hbmin, min_comp, consensus_name = "consensus"):
        def all_sequences_have_consensus_assigned(poagraph):
            return all([*map(lambda sequence: not sequence.active, poagraph.sources)])

        iteration_id = 0
        while not all_sequences_have_consensus_assigned(poagraph):
            try:
                new_poagraph = self._run_single_consensus_generation(consensus_output_dir, poagraph, hbmin, "consensus_" + str(iteration_id))
            except NoConsensusFound:
                 print("NO MORE CONSENSUSES FOUND")
                 break

            maximally_consensus_compatible_sources_IDs = self._get_compatible(new_poagraph.sources,
                                                                                new_poagraph.consensuses[0],
                                                                                0)
            sources_ID_map, nodes_ID_map = new_poagraph.deactivate_different_then(maximally_consensus_compatible_sources_IDs)

            try:
                narrowed_poagraph = self._run_single_consensus_generation(consensus_output_dir, new_poagraph, "consensus_narrowed_" + str(iteration_id))
            except NoConsensusFound:
                print("NO CONSENSUS FOR NARROWED POAGRAPH FOUND")
                break

            enhanced_consensus = narrowed_poagraph.consensuses[0]
            enhanced_consenssus_nodes_IDs = [nodes_ID_map[node_temp_ID] for node_temp_ID in enhanced_consensus.nodes_IDs]

            new_consensus_ID = len(poagraph.consensuses)
            poagraph.add_consensus(Consensus(ID=new_consensus_ID,
                                                 name=enhanced_consensus.name,
                                                 title=enhanced_consensus.title,
                                                 nodes_IDs=enhanced_consenssus_nodes_IDs))
            poagraph.calculate_compatibility_to_consensuses()

            good_consensus_compatible_sources_IDs = self._get_compatible(poagraph.sources,
                                                                            poagraph.consensuses[new_consensus_ID],
                                                                            min_comp)
            for source_ID, source in enumerate(poagraph.sources):
                if source_ID in good_consensus_compatible_sources_IDs:
                    poagraph.sources[source_ID].consensusID = new_consensus_ID


            poagraph.activate_sources_with_consensus_unassigned()
            iteration_id += 1


        return poagraph

    def _get_compatible(self, sources, consensus, min_comp):
        compatibilities = consensus.compatibility_to_sources
        max_compatibility = max(compatibilities)
        return [sourceID for sourceID, compatibility in enumerate(compatibilities) if
                max_compatibility-compatibility <= min_comp and sources[sourceID].active == True]

    def generate_visualization(self, consensuses_comparison=False, graph_visualization=False):
        print('Generate visualization...')
        for p in self.poagraphs:
            vizualization_output_dir = t.create_child_dir(p.path, "visualization")
            visualizator = POAGraphVisualizator(p, vizualization_output_dir, self.data_type)
            visualizator.generate(consensuses_comparison, graph_visualization)

class NoConsensusFound(Exception):
    pass