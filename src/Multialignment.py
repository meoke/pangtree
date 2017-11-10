from subprocess import run
import toolkit as t
import maf_reader as maf_reader
import po_reader as po_reader
from POAGraphVisualizator import POAGraphVisualizator
from Sequence import Consensus

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

    def generate_consensus(self, option, hbmin, min_comp, comp_range, tresholds):
        for i, p in enumerate(self.poagraphs):
            consensus_output_dir = t.create_child_dir(p.path, "consensus")
            if option is 0:
                print('Generate consensuses (hbmin=', str(hbmin) + ') in one iteration...')
                new_poagraph = self._run_single_consensus_generation(consensus_output_dir, p, hbmin)
                new_poagraph.path = p.path
                self.poagraphs[i] = new_poagraph
            elif option is 1:
                print('Generate consensuses (hbmin=', str(hbmin), ', min_comp=', str(min_comp), ') iteratively...')
                self._run_iterative_consensus_generation(consensus_output_dir, p, hbmin, min_comp)
            elif option is 2:
                print('Generate tree based consensus (hbmin=', str(hbmin), ', min_comp=', str(min_comp), ', range=', comp_range, ', tresholds=', tresholds)
                self._run_tree_consensus_generation(consensus_output_dir, p, hbmin, min_comp, comp_range, tresholds)


    def _run_single_consensus_generation(self, consensus_output_dir, poagraph, hbmin, consensus_name = "consensus"):
        print("PO generation")
        poagraph_as_po = poagraph.generate_po()

        file_name = t.join_path(consensus_output_dir, consensus_name)
        hb_file_name = t.change_file_extension(file_name,  '.hb')
        with open(file_name, 'w') as output_po_file:
            output_po_file.write(poagraph_as_po)
        print("Run poa")
        run(['../bin/poa', '-read_msa', file_name, '-hb', '-po', hb_file_name, '../bin/blosum80.mat', '-hbmin',
             str(hbmin)])

        print("Parse po to poagraph")
        new_poagraph = po_reader.parse_to_poagraph(hb_file_name, poagraph.path)

        print("Check if consensus found")
        if all([src.consensusID == -1 for src in new_poagraph.sources]):
            raise NoConsensusFound()

        print("Calculate compatibility")
        new_poagraph.calculate_compatibility_to_consensuses()
        return new_poagraph

    def _run_iterative_consensus_generation(self, consensus_output_dir, poagraph, hbmin, min_comp, consensus_name = "consensus"):
        def all_sequences_have_consensus_assigned(poagraph):
            return all([*map(lambda sequence: not sequence.active, poagraph.sources)])

        iteration_id = 0
        while not all_sequences_have_consensus_assigned(poagraph):
            min_comp *= 2
            try:
                print("Iteration ", str(iteration_id))
                new_poagraph = self._run_single_consensus_generation(consensus_output_dir=consensus_output_dir,
                                                                     poagraph=poagraph,
                                                                     hbmin=hbmin,
                                                                     consensus_name="consensus_" + str(iteration_id))
            except NoConsensusFound:
                 print("NO MORE CONSENSUSES FOUND")
                 max_compatibility = 0
                 best_consensus = -1
                 for i, source in enumerate(poagraph.sources):
                     if source.consensusID == -1:
                         for cons in poagraph.consensuses:
                             if cons.compatibility_to_sources[i] > max_compatibility:
                                 max_compatibility = cons.compatibility_to_sources[i]
                                 best_consensus = cons.currentID
                         poagraph.sources[i].consensusID = best_consensus
                         max_compatibility = 0
                         best_consensus = -1
                 break

            print("Get compatible for new consensus")
            maximally_consensus_compatible_sources_IDs = self._get_compatible(sources=new_poagraph.sources,
                                                                                consensus=new_poagraph.consensuses[0],
                                                                                min_comp=0,
                                                                                consensuses=new_poagraph.consensuses)
            print("Deactivate not compatible")
            sources_ID_map, nodes_ID_map = poagraph.deactivate_different_then(maximally_consensus_compatible_sources_IDs)


            try:
                print("Generate consensus for narrowed graph")
                narrowed_poagraph = self._run_single_consensus_generation(consensus_output_dir=consensus_output_dir,
                                                                          # poagraph=new_poagraph,
                                                                          poagraph=poagraph,
                                                                          hbmin=hbmin,
                                                                          consensus_name="consensus_narrowed_" + str(iteration_id))
            except NoConsensusFound:
                print("NO CONSENSUS FOR NARROWED POAGRAPH FOUND")
                break

            enhanced_consensus = narrowed_poagraph.consensuses[0]
            enhanced_consenssus_nodes_IDs = [nodes_ID_map[node_temp_ID] for node_temp_ID in enhanced_consensus.nodes_IDs]

            new_consensus_ID = len(poagraph.consensuses)
            poagraph.add_consensus(Consensus(currentID=new_consensus_ID,
                                                 name="CONSENS"+str(new_consensus_ID),
                                                 title=enhanced_consensus.title,
                                                 nodes_IDs=enhanced_consenssus_nodes_IDs))
            print("Calculate new compatibility to consensuses")
            poagraph.calculate_compatibility_to_consensuses()

            print("Get compatible to new consensuses")
            good_consensus_compatible_sources_IDs = self._get_compatible(poagraph.sources,
                                                                            poagraph.consensuses[new_consensus_ID],
                                                                            min_comp,
                                                                            poagraph.consensuses)
            if not good_consensus_compatible_sources_IDs:
                print("Nothing compatible!!!")
                break

            print("Assign consensus to compatible")
            for source_ID, source in enumerate(poagraph.sources):
                if source_ID in good_consensus_compatible_sources_IDs:
                    poagraph.sources[source_ID].consensusID = new_consensus_ID

            print("Activate sources with consensus unassigned")
            poagraph.activate_sources_with_consensus_unassigned()
            iteration_id += 1


        return poagraph

    def _run_tree_consensus_generation(consensus_output_dir, p, hbmin, min_comp, comp_range, tresholds):
        print("Tree consensus...")

    def _get_compatible(self, sources, consensus, min_comp, consensuses):
        def mean(numbers):
            return float(sum(numbers)) / max(len(numbers), 1)

        def is_best_compatibility_for_source(sourceID, current_compatibility):
            for consensus in consensuses:
                if consensus.compatibility_to_sources[sourceID] > current_compatibility:
                    return False
            return True

        compatibilities = consensus.compatibility_to_sources
        max_compatibility = max(compatibilities)
        mean_compatibility = mean(compatibilities)
        return [sourceID for sourceID, compatibility in enumerate(compatibilities) if
                abs(max_compatibility-compatibility) <= mean_compatibility*min_comp and
                is_best_compatibility_for_source(sourceID, compatibility)]

    def generate_visualization(self, consensuses_comparison=False, graph_visualization=False):
        print('Generate visualization...')
        for p in self.poagraphs:
            vizualization_output_dir = t.create_child_dir(p.path, "visualization")
            visualizator = POAGraphVisualizator(p, vizualization_output_dir, self.data_type)
            visualizator.generate(consensuses_comparison, graph_visualization)

class NoConsensusFound(Exception):
    pass