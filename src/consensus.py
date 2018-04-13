from subprocess import run
import numpy as np
from Errors import NoConsensusFound, NoTresholdFound, StopExceeded
from POAGraphRef import POAGraphRef
import toolkit as t
import po_reader as po_reader
import po_writer as po_writer
import time

class SourceCompatibilitesOptions(object):
    def __init__(self, assigned_consensus_ID, possible_consensusID_to_compatibility):
        assigned_consensusID = assigned_consensus_ID
        possible_consensusID_to_compatibility = possible_consensusID_to_compatibility

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

    src_to_nodeID = {}
    smallest_comp_up_to_now = 10
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
        cutoff_for_node = _find_cutoff_new(comp_to_current_srcs,
                                                multiplier,
                                                smallest_comp_up_to_now)

        #7 get sources compatible enough to the top consensus
        compatible_sources_IDs = get_compatible(current_srcs, comp_to_current_srcs, cutoff_for_node)

        #8 check if splitting should be continued based on stop condition - moved to the top o

        # wersja z próbą utrzymania węzłów na podobnym poziomie
        #TODOCZWARTEK2
            # decide how many sources will be added to the consensus
            # if children_nodes_IDs:
            #     the_smallest_comp_up_to_now = poagraph.get_min_cutoff(children_nodes_IDs) #todo a może jednak max
            # else:
            #     the_smallest_comp_up_to_now = 1

            # if the_smallest_comp_up_to_now < cutoff_for_node:
            #     max_c, max_consensus_nodes = get_top_consensus(poagraph, current_srcs, hbmin)
            #     konsensus_comps = [poagraph.get_comp(max_consensus_nodes, src_ID) for src_ID in current_srcs]
            #     srcs_to_include = current_srcs
            #     new_children_node_comp = min(konsensus_comps)
            #     current_srcs = []

                # srcs_to_include = current_srcs
                # current_srcs = []
                # new_children_node_comp = min(comp_to_current_srcs)
            # else:
        #wersja bez tego (w razie wrócenia do zakomentowanej - wsunąć do elsa

        srcs_to_include = compatible_sources_IDs
        current_srcs = np.setdiff1d(current_srcs, compatible_sources_IDs)
        new_children_node_comp = cutoff_for_node
        if new_children_node_comp < smallest_comp_up_to_now:
            smallest_comp_up_to_now = new_children_node_comp
        #### koniec przesuniecia


        # add the consensus
        poagraph.add_consensus(max_c, max_consensus_nodes)

        new_node = POAGraphRef(parent_ID=tree_node_ID,
                               sources_IDs=srcs_to_include,
                               consensus_ID=poagraph.consensuses[-1].ID,
                               min_compatibility=new_children_node_comp)

        new_node_ID = poagraph.add_poagraphref(new_node, tree_node_ID)
        children_nodes_IDs.append(new_node_ID)

        for srcID in srcs_to_include:
            src_to_nodeID[srcID] = new_node_ID

    if re_consensus:
        pass
        #then check if some sequences should be reassigned

    #TODOCZWARTEK1
    # for srcID, assigned_consensusID in src_to_nodeID.items():
    #     node_ID_to_comp = {}
    #     for child_node_ID in children_nodes_IDs:
    #         child_node = poagraph._poagraphrefs[child_node_ID]
    #         consensus_ID = child_node.consensus_ID
    #         possible_src_comp = poagraph.consensuses[consensus_ID].compatibility_to_sources[srcID]
    #         node_ID_to_comp[child_node_ID] = possible_src_comp
    #
    #     current_src_comp = poagraph.consensuses[assigned_consensusID].compatibility_to_sources[srcID]
    #     the_best_comp = max(node_ID_to_comp.values())
    #     # (the_best_comp, the_best_nodeID) =
    #     # if the_best_comp > current_src_comp:
    #         # usun z obecnego węzła
    #         # dodaj do nowego węzła
    #
    #
    #
    # for srcID in tree_node_src_IDs:
    #     node_ID_to_comp = {}
    #     # check consensus ID for the src
    #     for child_node_ID in children_nodes_IDs:
    #         child_node = poagraph[child_node_ID]
    #         consensus_ID = child_node.consensus_ID
    #         if srcID in child_node.sources_IDs:
    #             current_src_comp = poagraph.consensuses[consensus_ID].compatibility_to_sources[srcID]
    #         else:
    #             possible_src_comp = poagraph.consensuses[consensus_ID].compatibility_to_sources[srcID]
    #             node_ID_to_comp[child_node_ID] = possible_src_comp
    #
    #     # check if there is better comp in children nodes ids:
    #     the_best_comp = max(node_ID_to_comp.values())
    #     if current_src_comp < the_best_comp:
    #         for nodeID, comp in node_ID_to_comp.items():
    #             if comp == the_best_comp:
    #                 better_node_ID = nodeID
    #         poagraph._poagraphrefs[better_node_ID].sources_IDs.append(srcID)
            # poagraph._poagraphrefs[]
        # consensus_ID =
        # src_comp_to_currently_assigned_node =

    return children_nodes_IDs


def get_top_consensus(poagraph, sources_IDs, hbmin):
    po_file_path, nodes_map = po_writer.save_as_po(poagraph, sources_IDs)
    hb_file_path = t.change_file_extension(po_file_path, '.hb')
    run(['../bin/poa', '-read_msa', po_file_path, '-hb', '-po', hb_file_path, '../bin/blosum80.mat', '-v', '-hbmin',
         str(hbmin)])
    try:
        consensus0, consensus_nodes = po_reader.read_consensus(hb_file_path, consensusID=0)
    except NoConsensusFound:
        raise NoConsensusFound()

    consensus_actual_nodes = np.zeros(shape=len(poagraph.nodes), dtype=np.bool)
    for i, val in enumerate(consensus_nodes):
        orig_ID = nodes_map['orig_ID'][nodes_map['temp_ID'] == i][0]
        consensus_actual_nodes[orig_ID] = val

    return consensus0, consensus_actual_nodes


def get_compatible(poagraphref_srcs_IDs, compatibilities, cutoff_for_node):
    #todo użyć re_consensensusu
    indexes = np.where(compatibilities>=cutoff_for_node)
    return poagraphref_srcs_IDs[indexes]

def _find_cutoff_new(compatibilities, multiplier, smallest_comp_up_to_now):
    #dołożenie najmniejszego comp do tej pory
    values_to_search_for_cutoff = compatibilities.copy()
    if smallest_comp_up_to_now != 10:
        values_to_search_for_cutoff.append(smallest_comp_up_to_now)

    #policzyć średnie odległości na całej osi
    sorted_compatibilities = sorted(values_to_search_for_cutoff)
    distances = [abs(sorted_compatibilities[i + 1] - sorted_compatibilities[i]) for i in range(len(sorted_compatibilities) - 1)]
    mean_distance = t.mean(distances)

    #sprawdzić pierwsze przekroczenie na kawałku [0;smallest_comp_up_to_now]
    if smallest_comp_up_to_now != 10:
        range_end = sorted_compatibilities.index(smallest_comp_up_to_now)
    else:
        range_end = len(sorted_compatibilities)

    level_boundary = mean_distance * multiplier
    for i in range(len(sorted_compatibilities[0:range_end]) - 1):
        if sorted_compatibilities[i + 1] - sorted_compatibilities[i] >= level_boundary:
                return sorted_compatibilities[i + 1]

    if range_end < len(sorted_compatibilities):
        return sorted_compatibilities[range_end]
    else:
        return sorted_compatibilities[-1]


        #if nie ma:
        # return pierwszą wartość za smallest comp
        #else:
        # return prawą granicę przedziału z przekroczeniem


def _find_cutoff_for_node(compatibilities, multiplier, smallest_comp_up_to_now):
    #todo przyjrzeć się i zrobić testy
    sorted_compatibilities = sorted(compatibilities)
    distances = [abs(compatibilities[i+1] - compatibilities[i]) for i in range(len(compatibilities)-1)]
    if not distances:
        return compatibilities[0]
    mean_distance = t.mean(distances)
    level_boundary = mean_distance * multiplier

    for i in range(len(compatibilities)-1):
        if sorted_compatibilities[i+1] - sorted_compatibilities[i] >= level_boundary:
            return sorted_compatibilities[i+1]
    print("NO TRESHOLD FOUND UWAGA")
    print("Compatibilities:")
    print(compatibilities)
    print("Multiplier:")
    print(multiplier)
    print("Returned sorted_compatibilities[0]:")
    print(sorted_compatibilities[0])
    return sorted_compatibilities[0]
    raise NoTresholdFound()


def _find_max_cutoff(compatibilities, cutoff_search_range):
    #todo przeanalizować dokładniej i zrobić testy
    #todo a może wyrzucić cutoff search range i po prostu szukać największej zmiany? porównać wyniki.
    min_search_idx = round((len(compatibilities)-1) * cutoff_search_range[0]/100)
    max_search_idx = round((len(compatibilities)-1) * cutoff_search_range[1]/100)
    compatibilities_to_be_searched = sorted(compatibilities)[min_search_idx: max_search_idx]

    max_diff = 0
    if min_search_idx == max_search_idx == len(compatibilities)-1:
        return sorted(compatibilities)[-1]
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