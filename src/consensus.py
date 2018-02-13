from subprocess import run
import numpy as np
from Errors import NoConsensusFound, NoTresholdFound, StopExceeded
from POAGraphRef import POAGraphRef
import toolkit as t
import po_reader as po_reader
import po_writer as po_writer


def process_tree_node(poagraph, tree_node_ID, cutoff_search_range, multiplier, re_consensus, stop):
    def cancel_splitting(message):
        raise StopExceeded(message)

    hbmin = 0.2
    tree_node_compatibility = poagraph.get_poagraphref_compatibility(tree_node_ID)
    tree_node_src_IDs = poagraph.get_poagraphref_sources_IDs(tree_node_ID)
    if tree_node_compatibility >= stop or len(tree_node_src_IDs) == 1:
        return []

    current_srcs = tree_node_src_IDs
    children_nodes_IDs = []
    while len(current_srcs):
        #1 find top consensus for all sources
        c, c_nodes = get_top_consensus(poagraph, current_srcs, hbmin)

        #2 get compatibility of the top consensus to all sources in currently processed sources
        comp_to_current_srcs = [poagraph.get_comp(c_nodes, src_ID) for src_ID in current_srcs]

        #3 get cutoff based on search range
        max_cutoff = _find_max_cutoff(comp_to_current_srcs, cutoff_search_range)

        #4 get sources maximally compatible to top consensus and find consensus for them
        max_compatible_sources_IDs = current_srcs[np.where(comp_to_current_srcs>=max_cutoff)]
        max_c, max_consensus_nodes = get_top_consensus(poagraph, max_compatible_sources_IDs, hbmin)

        #5 get compatibility of max consensus to all current sources
        comp_to_current_srcs = [poagraph.get_comp(max_consensus_nodes, src_ID) for src_ID in current_srcs]

        #6 get cutoff based on compatibilties
        cutoff_for_node = _find_cutoff_for_node(comp_to_current_srcs, multiplier)

        #7 get sources compatible enough to the top consensus
        compatible_sources_IDs = get_compatible(current_srcs, comp_to_current_srcs, cutoff_for_node,
                                                re_consensus)

        #8 check if splitting should be continued based on stop condition - moved to the top
        # if cutoff_for_node > stop:
        #     cancel_splitting("Stop was exceeded.")

        # add the consensus
        poagraph.add_consensus(max_c, max_consensus_nodes)

        # decide how many sources will be added to the consensus
        if children_nodes_IDs:
            the_smallest_comp_up_to_now = poagraph.get_min_cutoff(children_nodes_IDs) #todo a może jednak max
        else:
            the_smallest_comp_up_to_now = 1

        if the_smallest_comp_up_to_now < cutoff_for_node:
            srcs_to_include = current_srcs
            current_srcs = []
            new_children_node_comp = min(comp_to_current_srcs)
        else:
            srcs_to_include = compatible_sources_IDs
            current_srcs = np.setdiff1d(current_srcs, compatible_sources_IDs)
            new_children_node_comp = cutoff_for_node

        new_node = POAGraphRef(parent_ID=tree_node_ID,
                               sources_IDs=srcs_to_include,
                               consensus_ID=poagraph.consensuses[-1].ID,
                               min_compatibility=new_children_node_comp)
        new_node_ID = poagraph.add_poagraphref(new_node, tree_node_ID)
        children_nodes_IDs.append(new_node_ID)

    return children_nodes_IDs


def get_top_consensus(poagraph, sources_IDs, hbmin):
    po_file_path, nodes_map = po_writer.save_as_po(poagraph, sources_IDs)
    hb_file_path = t.change_file_extension(po_file_path, '.hb')
    run(['../bin/poa', '-read_msa', po_file_path, '-hb', '-po', hb_file_path, '../bin/blosum80.mat', '-hbmin',
         str(hbmin)])
    try:
        consensus0, consensus_nodes = po_reader.read_consensus(hb_file_path, consensusID=0)
    except NoConsensusFound:
        raise NoConsensusFound()

    consensus_actual_nodes = np.zeros(shape=len(poagraph.nodes), dtype=np.bool)
    node_ID = 0
    for i, val in enumerate(consensus_nodes):
        orig_ID = nodes_map['orig_ID'][nodes_map['temp_ID'] == i][0]
        consensus_actual_nodes[orig_ID] = val
        #node_ID = node_ID + 1

    return consensus0, consensus_actual_nodes


def get_compatible(poagraphref_srcs_IDs, compatibilities, cutoff_for_node, re_consensus):
    return poagraphref_srcs_IDs[np.where(compatibilities>=cutoff_for_node)]

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


def _find_cutoff_for_node(compatibilities, multiplier):
    #todo przyjrzeć się
    sorted_compatibilities = sorted(compatibilities)
    distances = [abs(compatibilities[i+1] - compatibilities[i]) for i in range(len(compatibilities)-1)]
    if not distances:
        return compatibilities[0]
    mean_distance = t.mean(distances)
    level_boundary = mean_distance * multiplier

    for i in range(len(compatibilities)-1):
        if sorted_compatibilities[i+1] - sorted_compatibilities[i] >= level_boundary:
            return sorted_compatibilities[i+1]
    raise NoTresholdFound()


def _find_max_cutoff(compatibilities, cutoff_search_range):
    #todo przeanalizować dokładniej
    min_search_idx = round((len(compatibilities)-1) * cutoff_search_range[0]/100)
    max_search_idx = round((len(compatibilities)-1) * cutoff_search_range[1]/100)
    compatibilities_to_be_searched = sorted(compatibilities)[min_search_idx: max_search_idx]

    max_diff = 0
    if min_search_idx == max_search_idx == len(compatibilities)-1:
        return sorted(compatibilities)[-1]#sorted(compatibilities)[min_search_idx-1]
    elif min_search_idx == max_search_idx:
        return sorted(compatibilities)[min_search_idx]
    elif not compatibilities_to_be_searched:
        raise ValueError("No compatibilites to be searched.")
    else:
        cutoff_value = compatibilities_to_be_searched[0]

    for i, comp in enumerate(compatibilities_to_be_searched):
        if i < (len(compatibilities_to_be_searched) - 1) and compatibilities_to_be_searched[i + 1] - comp > max_diff:
            max_diff = compatibilities_to_be_searched[i + 1] - comp
            cutoff_value = compatibilities_to_be_searched[i + 1]
    return cutoff_value