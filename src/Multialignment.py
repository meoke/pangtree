from subprocess import run
import toolkit as t
import maf_reader as maf_reader
import po_reader as po_reader
#from POAGraphVisualizator import POAGraphVisualizator
#from Sequence import Consensus
#from Errors NoConsensusFound import *
#from NoTresholdFound import *
from POAGraphRef import POAGraphRef
import numpy as np
import consensus as cons


class Multialignment(object):
    def __init__(self, data_type='ebola'):
        self.name = None
        self.output_dir = None
        self.poagraphs = None
        self.data_type = data_type
        #self.tresholds = [] #TODO to jest wlasiciwie niepotrzebne, skoro poagraphrefs ma to ogarniac

    def build_multialignment_from_maf(self, maf_file_path, merge_option):
        """Build multialignment structrure from maf file.

        Keyword arguments:
        maf_file_path -- path to the maf file
        merge_option -- #todo uzupelnic poprawnie
        """
        print("Building multialignment from " + maf_file_path + "...") #todo logging

        self.name = self._get_multialignment_name(maf_file_path)
        self.output_dir = self._get_output_dir(maf_file_path)
        self.poagraphs = maf_reader.parse_to_poagraphs(file_path=maf_file_path,
                                                       merge_option=merge_option,
                                                       multialignment_name=self.name,
                                                       output_dir=self.output_dir)

    def build_multialignment_from_po(self, po_file_name):
        print("Building multialignment from " + po_file_name + "...") #todo logging

        self.name = self._get_multialignment_name(po_file_name)
        self.output_dir = self._get_output_dir(po_file_name)
        self.poagraphs = [po_reader.parse_to_poagraph(file_path=po_file_name,
                                                      output_dir=self.output_dir)]

    def _get_multialignment_name(self, input_file_name):
        return t.get_file_name_without_extension(input_file_name)

    def _get_output_dir(self, input_file_name):
        return t.create_next_sibling_dir(input_file_name, "converted")

    def generate_consensus(self, option, hbmin, min_comp, comp_range, tresholds, multiplier, stop, re_consensus):
        for i, p in enumerate(self.poagraphs):
            #consensus_output_dir = t.create_child_dir(p.path, "consensus")
            #visualization_output_dir = t.create_child_dir(p.path, "visualization")
            if option is 1:
                print('Generate consensuses (hbmin=', str(hbmin) + ') in one iteration...') # todo logging
                # new_poagraph = self._run_single_consensus_generation(consensus_output_dir, p, hbmin)
                # new_poagraph.path = p.path
                # self.poagraphs[i] = new_poagraph
            # elif option is 2:
            #     print('Generate consensuses (hbmin=', str(hbmin), ', min_comp=', str(min_comp), ') iteratively...')
            #     self._run_iterative_consensus_generation(consensus_output_dir, p, hbmin, min_comp)
            elif option is 3:
                #print('Generate tree based consensus (hbmin=', str(hbmin), ', min_comp=', str(min_comp), ', range=', comp_range, ', tresholds=', tresholds)
                cutoff_search_range = self._convert_str_to_tuple(comp_range)
                #tresholds = self._convert_str_to_list(tresholds)
                #self._run_tree_consensus_generation(consensus_output_dir, visualization_output_dir, p, hbmin, min_comp, cutoff_search_range, multiplier, stop, re_consensus)
                for poagraph in self.poagraphs:
                    consensus_output_dir = t.create_child_dir(poagraph.path, "tconsensus")
                    hbmin = 0.2

                    parent_tree_node = POAGraphRef(sources_IDs=np.array(range(len(poagraph.sources))))
                    tree_nodes_to_process = [parent_tree_node]

                    finished_sources_IDs = np.empty(len(poagraph.sources), dtype=np.int16)
                    finished_sources_IDs.fill(-1)
                    while -1 in finished_sources_IDs:
                        new_tree_nodes_to_process = []
                        for tree_node in tree_nodes_to_process: #poagraphref, hbmin, comp_range, multiplier, output_dir):
                            new_tree_nodes_to_process = cons.process_tree_node(poagraph,
                                                                               tree_node,
                                                                               hbmin,
                                                                               cutoff_search_range,
                                                                               multiplier,
                                                                               consensus_output_dir)
    #
    # def _run_single_consensus_generation(self, consensus_output_dir, poagraph, hbmin, consensus_name = "consensus"):
    #     print('PO generation')
    #     poagraph_as_po = poagraph.generate_po()
    #
    #     file_name = t.join_path(consensus_output_dir, consensus_name)
    #     hb_file_name = t.change_file_extension(file_name,  '.hb')
    #     with open(file_name, 'w') as output_po_file:
    #         output_po_file.write(poagraph_as_po)
    #     print('Run poa')
    #     run(['../bin/poa', '-read_msa', file_name, '-hb', '-po', hb_file_name, '../bin/blosum80.mat', '-hbmin',
    #          str(hbmin)])
    #
    #     print('Parse po to poagraph')
    #     new_poagraph = po_reader.parse_to_poagraph(hb_file_name, poagraph.path)
    #
    #     print('Check if consensus found')
    #     if all([src.consensusID == -1 for src in new_poagraph.sources]):
    #         raise NoConsensusFound()
    #
    #     print('Calculate compatibility')
    #     for source in new_poagraph.sources:
    #         source.consensuses = [source.consensusID]
    #
    #     new_poagraph.calculate_compatibility_to_consensuses()
    #     return new_poagraph
    #
    # def _run_iterative_consensus_generation(self, consensus_output_dir, poagraph, hbmin, min_comp, consensus_name = "consensus"):
    #     # todo naprawic pod katem source.consensuses
    #     # todo czy to w ogole zostawiac?
    #     def all_sequences_have_consensus_assigned(poagraph):
    #         return all([*map(lambda sequence: not sequence.active, poagraph.sources)])
    #
    #     iteration_id = 0
    #     while not all_sequences_have_consensus_assigned(poagraph):
    #         min_comp *= 2
    #         try:
    #             print("Iteration ", str(iteration_id))
    #             new_poagraph = self._run_single_consensus_generation(consensus_output_dir=consensus_output_dir,
    #                                                                  poagraph=poagraph,
    #                                                                  hbmin=hbmin,
    #                                                                  consensus_name="consensus_" + str(iteration_id))
    #         except NoConsensusFound:
    #              print("NO MORE CONSENSUSES FOUND")
    #              max_compatibility = 0
    #              best_consensus = -1
    #              for i, source in enumerate(poagraph.sources):
    #                  if source.consensusID == -1:
    #                      for cons in poagraph.consensuses:
    #                          if cons.compatibility_to_sources[i] > max_compatibility:
    #                              max_compatibility = cons.compatibility_to_sources[i]
    #                              best_consensus = cons.currentID
    #                      poagraph.sources[i].consensusID = best_consensus
    #                      max_compatibility = 0
    #                      best_consensus = -1
    #              break
    #
    #         print("Get compatible for new consensus")
    #         maximally_consensus_compatible_sources_IDs = self._get_compatible(sources=new_poagraph.sources,
    #                                                                             consensus=new_poagraph.consensuses[0],
    #                                                                             min_comp=0,
    #                                                                             consensuses=new_poagraph.consensuses)
    #         print("Deactivate not compatible")
    #         sources_ID_map, nodes_ID_map = poagraph.deactivate_different_then(maximally_consensus_compatible_sources_IDs)
    #
    #
    #         try:
    #             print("Generate consensus for narrowed graph")
    #             narrowed_poagraph = self._run_single_consensus_generation(consensus_output_dir=consensus_output_dir,
    #                                                                       # poagraph=new_poagraph,
    #                                                                       poagraph=poagraph,
    #                                                                       hbmin=hbmin,
    #                                                                       consensus_name="consensus_narrowed_" + str(iteration_id))
    #         except NoConsensusFound:
    #             print("NO CONSENSUS FOR NARROWED POAGRAPH FOUND")
    #             break
    #
    #         enhanced_consensus = narrowed_poagraph.consensuses[0]
    #         enhanced_consenssus_nodes_IDs = [nodes_ID_map[node_temp_ID] for node_temp_ID in enhanced_consensus.nodes_IDs]
    #
    #         new_consensus_ID = len(poagraph.consensuses)
    #         poagraph.add_consensus(Consensus(currentID=new_consensus_ID,
    #                                              name="CONSENS"+str(new_consensus_ID),
    #                                              title=enhanced_consensus.title,
    #                                              nodes_IDs=enhanced_consenssus_nodes_IDs))
    #         print("Calculate new compatibility to consensuses")
    #         poagraph.calculate_compatibility_to_consensuses()
    #
    #         print("Get compatible to new consensuses")
    #         good_consensus_compatible_sources_IDs = self._get_compatible(poagraph.sources,
    #                                                                         poagraph.consensuses[new_consensus_ID],
    #                                                                         min_comp,
    #                                                                         poagraph.consensuses)
    #         if not good_consensus_compatible_sources_IDs:
    #             print("Nothing compatible!!!")
    #             break
    #
    #         print("Assign consensus to compatible")
    #         for source_ID, source in enumerate(poagraph.sources):
    #             if source_ID in good_consensus_compatible_sources_IDs:
    #                 poagraph.sources[source_ID].consensusID = new_consensus_ID
    #
    #         print("Activate sources with consensus unassigned")
    #         poagraph.activate_sources_with_consensus_unassigned()
    #         iteration_id += 1
    #
    #
    #     return poagraph
    #
    # def _run_tree_consensus_generation(self, consensus_output_dir, visualization_output_dir, poagraph, hbmin, min_comp, comp_range, multiplier, stop, re_consensus):
    #     def current_comp_higher_then_earlier(comp, srcID, other_poagraphRefs):
    #         return True
    #         # for poagraphRef in other_poagraphRefs:
    #         #     if poagraphRef.consensus.compatibility_to_sources[srcID] > comp:
    #         #         return False
    #         #     else:
    #         #         try:
    #         #             poagraphRef.sourcesIDs.remove(srcID)
    #         #         except ValueError:
    #         #             continue
    #         # return True
    #
    #     def get_not_yet_classified_sources_IDs(poagraphRef, current_treshold_poagraphRefs):
    #         classified_sources = [srcID for poagraphRef in current_treshold_poagraphRefs for srcID in poagraphRef.sourcesIDs]
    #         all_sources = [sourceID for sourceID in poagraphRef.sourcesIDs]
    #         return sorted(list(set(all_sources) - set(classified_sources)))
    #
    #     def not_all_sources_finished(finished_sources):
    #         if sorted(finished_sources) != [src.currentID for src in poagraph.sources]:
    #             return True
    #         return False
    #
    #     def find_parent_consensus_id(poagraph, example_source_ID):
    #         for consensus in reversed(poagraph.consensuses):
    #             if example_source_ID in consensus.sources_IDs:
    #                 return consensus.currentID
    #         return -1
    #
    #     hbmin = 0.2
    #     poagraphRefs = [POAGraphRef([src.currentID for src in poagraph.sources])]
    #     current_treshold_poagraphRefs = []
    #     finished_sources = []
    #     #for tr, treshold in enumerate(tresholds):
    #     while not_all_sources_finished(finished_sources):
    #         print("TRESHOLD NUMBER " + str(len(self.tresholds)))
    #
    #         for pr, poagraphRef in enumerate(poagraphRefs):
    #             iteration_id = -1
    #             sourcesIDs_to_classify = poagraphRef.sourcesIDs
    #             poagraphRef_tresholds = []
    #             while sourcesIDs_to_classify:
    #                 iteration_id += 1
    #
    #                 try:
    #                     consensus_from_all = self._get_top_consensus_for_poagraph_part(consensus_output_dir=consensus_output_dir,
    #                                                           poagraph=poagraph,
    #                                                           sourcesIDs_to_be_used=sourcesIDs_to_classify,
    #                                                           hbmin=hbmin,
    #                                                           consensus_name="all_consensus_" + str(len(self.tresholds)) + '_' + str(pr) + '_' + str(iteration_id))
    #                 except NoConsensusFound:
    #                      print("NO MORE CONSENSUSES FOUND")
    #                      break
    #
    #                 cons_from_all_comp_to_srcs_IDs_to_classify = [comp for srcID, comp in
    #                                                                                enumerate(consensus_from_all.compatibility_to_sources)
    #                                                                                if consensus_from_all.sources_IDs[srcID] in sourcesIDs_to_classify]
    #
    #                 cutoff_value = self._find_cutoff_value(cons_from_all_comp_to_srcs_IDs_to_classify, comp_range)
    #
    #
    #                 maximally_consensus_compatible_sources_IDs = [consensus_from_all.sources_IDs[srcID] for srcID, comp in
    #                                                               enumerate(consensus_from_all.compatibility_to_sources) if
    #                                                               comp >= cutoff_value and
    #                                                               consensus_from_all.sources_IDs[srcID] in sourcesIDs_to_classify]
    #
    #                 try:
    #                     top_consensus = self._get_top_consensus_for_poagraph_part(consensus_output_dir=consensus_output_dir,
    #                                                                               poagraph=poagraph,
    #                                                                               sourcesIDs_to_be_used=maximally_consensus_compatible_sources_IDs,
    #                                                                               hbmin=hbmin,
    #                                                                               consensus_name="top_consensus_" + str(len(self.tresholds)) + '_' + str(pr) + '_' + str(iteration_id), sourcesIDs_to_calc_compatbility = sourcesIDs_to_classify)
    #                 except NoConsensusFound:
    #                      print("NO TOP CONSENSUSES FOUND")
    #
    #                 if poagraphRef_tresholds and all(map(lambda comp : min(poagraphRef_tresholds) < comp, top_consensus.compatibility_to_sources)):
    #                     consensus_compatible_sources_IDs = sourcesIDs_to_classify
    #                     if min(top_consensus.compatibility_to_sources) >= stop:
    #                         finished_sources += consensus_compatible_sources_IDs
    #                         print("Treshold exceeded stop value.")
    #                 else:
    #                     treshold = self._find_current_treshold(top_consensus.compatibility_to_sources, multiplier)
    #                     poagraphRef_tresholds.append(treshold)
    #                     self.tresholds.append(min(top_consensus.compatibility_to_sources))
    #                     consensus_compatible_sources_IDs = [top_consensus.sources_IDs[srcID] for srcID, comp in enumerate(top_consensus.compatibility_to_sources) if
    #                                                     comp >= treshold and
    #                                                     current_comp_higher_then_earlier(comp, srcID, current_treshold_poagraphRefs)]
    #
    #                     if treshold >= stop:
    #                         finished_sources += consensus_compatible_sources_IDs
    #                         print("Treshold exceeded stop value.")
    #
    #                 top_consensus.compatibility_to_sources = [comp for i, comp in enumerate(top_consensus.compatibility_to_sources) if top_consensus.sources_IDs[i] in consensus_compatible_sources_IDs]
    #                 top_consensus.sources_IDs = consensus_compatible_sources_IDs
    #
    #                 #dołożenie żeby były wszystkie consensusy w tym poagrafie?
    #                 ##
    #
    #                 new_poagraphRef = POAGraphRef(consensus_compatible_sources_IDs, top_consensus)
    #                 current_treshold_poagraphRefs.append(new_poagraphRef)
    #                 sourcesIDs_to_classify = get_not_yet_classified_sources_IDs(poagraphRef, current_treshold_poagraphRefs)
    #
    #         # enumerate_consensuses(current_treshold_poagraphRefs)
    #         for i, poagraphRef in enumerate(current_treshold_poagraphRefs):
    #             consensus = poagraphRef.consensus
    #             consensusID = len(poagraph.consensuses)
    #             consensus.currentID = consensusID
    #             consensus.name = "CONSENS" + str(consensusID)
    #             # consensus.level = len(self.tresholds)
    #             consensus.parent_consensus = find_parent_consensus_id(poagraph, poagraphRef.sourcesIDs[0])
    #             consensus.level = min(consensus.compatibility_to_sources)
    #             consensus.sources_IDs = poagraphRef.sourcesIDs ##DODANE
    #             poagraph.add_consensus(consensus)
    #             poagraph.calculate_compatibility_to_consensuses(consensusID=len(poagraph.consensuses) - 1)
    #             for srcID in poagraphRef.sourcesIDs:
    #                 poagraph.sources[srcID].consensuses[consensus.level] = consensusID
    #         poagraphRefs = []
    #         for p in current_treshold_poagraphRefs:
    #             if p.sourcesIDs[0] in finished_sources:
    #                 continue
    #             else:
    #                 poagraphRefs.append(p)
    #         # poagraphRefs = current_treshold_poagraphRefs
    #         current_treshold_poagraphRefs = []
    #
    #     # final consensuses assignment
    #     for i, poagraphRef in enumerate(poagraphRefs):
    #         for srcID in poagraphRef.sourcesIDs:
    #             poagraph.sources[srcID].consensusID = i
    #
    # def _get_top_consensus_for_poagraph_part(self, consensus_output_dir, poagraph, sourcesIDs_to_be_used, hbmin, consensus_name="consensus", sourcesIDs_to_calc_compatbility = None):
    #     print("PO generation for particular sources IDs")
    #     poagraph_as_po, new_to_original_nodes_IDs = poagraph.generate_partial_po(sourcesIDs_to_use=sourcesIDs_to_be_used)
    #
    #     file_name = t.join_path(consensus_output_dir, consensus_name)
    #     hb_file_name = t.change_file_extension(file_name, '.hb')
    #     with open(file_name, 'w') as output_po_file:
    #         output_po_file.write(poagraph_as_po)
    #     print("Run poa")
    #     run(['../bin/poa', '-read_msa', file_name, '-hb', '-po', hb_file_name, '../bin/blosum80.mat',
    #          '-hbmin',
    #          str(hbmin)])
    #
    #     print("Parse po to poagraph")
    #     try:
    #         consensus0 = po_reader.read_single_consensus(hb_file_name, consensusID=0)
    #     except NoConsensusFound:
    #         raise NoConsensusFound()
    #
    #     for i, nodeID in enumerate(consensus0.nodes_IDs):
    #         consensus0.nodes_IDs[i] = new_to_original_nodes_IDs[consensus0.nodes_IDs[i]]
    #
    #     print("Calculate compatibility")
    #     poagraph.add_consensus(consensus0)
    #     if sourcesIDs_to_calc_compatbility:
    #         srcs = sourcesIDs_to_calc_compatbility
    #     else:
    #         srcs = sourcesIDs_to_be_used
    #     poagraph.calculate_compatibility_to_consensuses(consensusID=len(poagraph.consensuses)-1, sources_IDs=srcs)
    #     last_consensus = poagraph.consensuses[len(poagraph.consensuses)-1]
    #     poagraph.remove_last_consensus()
    #     last_consensus.sources_IDs = srcs
    #     return last_consensus
    #
    def _convert_str_to_tuple(self, comp_range):
        return eval(comp_range[1:-1])
    #
    # def _convert_str_to_list(self, tresholds):
    #     return eval(tresholds)
    #
    # def _find_cutoff_value(self, compatibilities, comp_range):
    #     min_idx = round((len(compatibilities)-1) * comp_range[0]/100)+1
    #     max_idx = round((len(compatibilities)-1) * comp_range[1]/100)+1
    #     compatibilities_to_be_searched = sorted(compatibilities)[min_idx: max_idx]
    #
    #     max_diff = 0
    #     if min_idx == max_idx == len(compatibilities):
    #         return sorted(compatibilities)[min_idx-1]
    #     elif min_idx == max_idx:
    #         return sorted(compatibilities)[min_idx]
    #     elif not compatibilities_to_be_searched:
    #         raise ValueError("No compatibilites to be searched.")
    #     else:
    #         cutoff_value = compatibilities_to_be_searched[0]
    #
    #     for i, comp in enumerate(compatibilities_to_be_searched):
    #         if i < (len(compatibilities_to_be_searched) - 1) and compatibilities_to_be_searched[i + 1] - comp > max_diff:
    #             max_diff = compatibilities_to_be_searched[i + 1] - comp
    #             cutoff_value = compatibilities_to_be_searched[i + 1]
    #     return cutoff_value
    #
    # def _get_compatible(self, sources, consensus, min_comp, consensuses):
    #     def mean(numbers):
    #         return float(sum(numbers)) / max(len(numbers), 1)
    #
    #     def is_best_compatibility_for_source(sourceID, current_compatibility):
    #         for consensus in consensuses:
    #             if consensus.compatibility_to_sources[sourceID] > current_compatibility:
    #                 return False
    #         return True
    #
    #     compatibilities = consensus.compatibility_to_sources
    #     max_compatibility = max(compatibilities)
    #     mean_compatibility = mean(compatibilities)
    #     return [sourceID for sourceID, compatibility in enumerate(compatibilities) if
    #             abs(max_compatibility-compatibility) <= mean_compatibility*min_comp and
    #             is_best_compatibility_for_source(sourceID, compatibility)]
    #
    # def _find_current_treshold(self, compatibilities, multiplier):
    #     sorted_compatibilities = sorted(compatibilities)
    #     distances = [abs(compatibilities[i+1] - compatibilities[i]) for i in range(len(compatibilities)-1)]
    #     if not distances:
    #         return compatibilities[0]
    #     mean_distance = t.mean(distances)
    #     level_boundary = mean_distance * multiplier
    #
    #     for i in range(len(compatibilities)-1):
    #         if sorted_compatibilities[i+1] - sorted_compatibilities[i] >= level_boundary:
    #             return sorted_compatibilities[i+1]
    #     raise NoTresholdFound()
    #
    # def generate_visualization(self, consensuses_comparison=False, graph_visualization=False, processing_time='', tresholds= '', consensus_algorithm=False):
    #     print('Generate visualization...')
    #     for p in self.poagraphs:
    #         vizualization_output_dir = t.create_child_dir(p.path, "visualization")
    #         visualizator = POAGraphVisualizator(p, vizualization_output_dir, self.data_type)
    #         # visualizator.generate(consensuses_comparison, graph_visualization, processing_time, self._convert_str_to_list(tresholds), consensus_algorithm)
    #         visualizator.generate(consensuses_comparison, graph_visualization, processing_time, self.tresholds, consensus_algorithm)
    #


