from subprocess import run
import numpy as np
from Errors import NoConsensusFound, NoTresholdFound
from POAGraphRef import POAGraphRef
import toolkit as t
import po_reader as po_reader
import po_writer as po_writer

def process_tree_node(poagraph, poagraphref, hbmin, cutoff_search_range, multiplier, output_dir, re_consensus, stop):
    children_nodes = []
    sources_IDs_to_process = poagraphref.sources_IDs
    #todo cos jest nie tak
    while len(sources_IDs_to_process):
        consensus, consensus_nodes = get_top_consensus(poagraph, sources_IDs_to_process, hbmin, output_dir)
        compatibilities_to_tree_node_sources = [poagraph.calc_compatibility(consensus_nodes, src_ID) for src_ID in sources_IDs_to_process]
        cutoff_for_max = _find_cutoff_for_max(compatibilities_to_tree_node_sources, cutoff_search_range)
        max_compatible_sources_IDs = sources_IDs_to_process[np.where(compatibilities_to_tree_node_sources>=cutoff_for_max)]
        max_consensus, max_consensus_nodes = get_top_consensus(poagraph, max_compatible_sources_IDs, hbmin, output_dir)
        compatibilities_to_tree_node_sources = [poagraph.calc_compatibility(max_consensus_nodes, src_ID) for src_ID in sources_IDs_to_process]

        cutoff_for_node = _find_cutoff_for_node(compatibilities_to_tree_node_sources, multiplier)
        compatible_sources_IDs = get_compatible(poagraphref, compatibilities_to_tree_node_sources, cutoff_for_node,
                                                re_consensus)
        if cutoff_for_node > stop:
            sources_IDs_to_process = np.setdiff1d(sources_IDs_to_process, compatible_sources_IDs)
            continue


        poagraph.add_consensus(max_consensus, max_consensus_nodes)
        #poagraph_ref = POAGraphRef(compatible_sources_IDs, poagraph.consensuses[-1].ID, min(compatibilities_to_tree_nodes_sources))

        # current_min_probability = min(compatibilities_to_tree_node_sources)
        the_smallest_compaibility_up_to_now = get_min_treshold(children_nodes) #children_nodes[-1].min_compatibility wwszystkie children nodes
        if the_smallest_compaibility_up_to_now < cutoff_for_node:
            poagraph_ref = POAGraphRef(sources_IDs_to_process, poagraph.consensuses[-1].ID, min(compatibilities_to_tree_node_sources))
            children_nodes.append(poagraph_ref)
            sources_IDs_to_process = []
        else: #dzielimy dalej
            poagraph_ref = POAGraphRef(compatible_sources_IDs, poagraph.consensuses[-1].ID, cutoff_for_node)
            children_nodes.append(poagraph_ref)
            sources_IDs_to_process = np.setdiff1d(sources_IDs_to_process, compatible_sources_IDs)
            # tree_node.sources_IDs = np.array([src_ID for src_ID in tree_node.sources_IDs if src_ID not in compatible_sources_IDs])

    return children_nodes


def get_top_consensus(poagraph, sources_IDs, hbmin, output_dir):
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


def calc_compatibility(self, consensus, tree_node):
    pass


def get_min_treshold(children_nodes):
    return min([p.min_compatibility for p in children_nodes]) if children_nodes else 1 #todo to trochę hack


def get_compatible(poagraphref, compatibilities, cutoff_for_node, re_consensus):
    return poagraphref.sources_IDs[np.where(compatibilities>=cutoff_for_node)]

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


def _find_cutoff_for_max(compatibilities, cutoff_search_range):
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