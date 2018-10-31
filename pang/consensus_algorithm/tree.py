from graph.Pangraph import Pangraph
from pathlib import Path
from copy import deepcopy
from metadata.MultialignmentMetadata import MultialignmentMetadata
from . import simple
import numpy as np
from consensus_algorithm.TreeConfig import TreeConfig
from .AlgorithmParams import AlgorithmParams
from consensus_data.SubPangraph import SubPangraph
from consensus_data.TreeConsensusManager import TreeConsensusManager
from consensus_data.ConsensusNode import ConsensusNode
from consensus_data.Errors import NoConsensus
from statistics import mean
from collections import deque

# import logging
# logger = logging.getLogger(__name__)

ap = AlgorithmParams()


def run(outputdir: Path, pangraph: Pangraph, config: TreeConfig, genomes_info: MultialignmentMetadata) -> Pangraph:
    ap.outputdir, ap.config, ap.genomes_info = outputdir, config, genomes_info

    # all_sequences_ids = pangraph.get_path_ids()

    # subpangraph: SubPangraph = SubPangraph(pangraph, all_sequences_ids)
    # cm = TreeConsensusManager(subpangraph.get_nodes_count())

    # root_node = ConsensusNode(sequences_names=list(all_sequences_ids))
    # cm.add_node(root_node, get_top_consensus(subpangraph))

    # root_consensusManager = produce_tree(subpangraph, cm) # todo z takimi parametrami, żeby był root z jednym consensusem
    root_consensusManager = produce_tree2(pangraph)  # todo z takimi parametrami, żeby był root z jednym consensusem

    pangraph.set_consensus_manager(root_consensusManager)
    return pangraph


def node_ready(node: ConsensusNode):
    min_own_comp = min(node.get_compatibilities_to_own_sources())
    if len(node.sequences_names) == 1:
        return True
    if min_own_comp >= ap.config.stop:
        return True
    return False


def produce_tree2(pangraph: Pangraph) -> TreeConsensusManager:
    all_sequences_names = pangraph.get_path_names()  #

    cm = TreeConsensusManager(max_nodes_count=pangraph.get_nodes_count())  # jego będę produkować
    root_node = ConsensusNode(sequences_names=list(all_sequences_names))
    root_pangraph = SubPangraph(pangraph, pangraph.get_path_names())
    cm.add_node(root_node, get_top_consensus(root_pangraph))

    nodes_to_process = deque([root_node])
    while nodes_to_process:
        subtree_root = nodes_to_process.pop()
        current_node_pangraph = SubPangraph(pangraph, subtree_root.sequences_names)
        children_nodes_manager = get_children_cm2(current_node_pangraph, subtree_root)  # TreeConsensusManager z węzłami-braćmi i odpowiadającymi im consensusami w wersji dla pangraohu ograniczonego do ścieżek w danym nodzie

        chidren_nodes = children_nodes_manager.get_nodes()
        if len(chidren_nodes) == 1:
            continue

        for child in chidren_nodes:
            # consensus_in_subpangraph = children_nodes.get_consensus(child.consensus_id)
            # consensus = current_node_pangraph.get_consensus_remapped_to_original_nodes(child.consensus_id) #przekazać consensus z managera lokalnego
            consensus = children_nodes_manager.get_consensus(child.consensus_id)
            child.parent_node_id = subtree_root.consensus_id
            child.compatibilities_to_all = pangraph.get_paths_compatibility_to_consensus(consensus)
            child_node_id = cm.add_node(child, consensus)
            subtree_root.children_nodes.append(child_node_id)

            if not node_ready(child):
                nodes_to_process.append(child)
    return cm

def get_children_cm2(subpangraph: SubPangraph, node: ConsensusNode) -> TreeConsensusManager:
    current_paths_names = node.sequences_names
    # current_paths_ids = [subpangraph.pangraph.get_path_id(path_name) for path_name in current_paths_names]
    # local_consensus_manager = TreeConsensusManager(max_nodes_count=subpangraph.get_nodes_count())
    local_consensus_manager = TreeConsensusManager(max_nodes_count=subpangraph.orig_nodes_count)
    s = deepcopy(subpangraph)
    while current_paths_names:
        subpangraph = run_poa(subpangraph)
        c_to_node = subpangraph.get_paths_compatibility(0) #zgodnie ze swoją przechowywaną kolejnością
        max_cutoff = find_max_cutoff(c_to_node, ap.config.cutoff_search_range)
        max_c_sources_names = get_max_compatible_sources_ids(current_paths_names, c_to_node, max_cutoff)

        subsubpangraph = SubPangraph(subpangraph.pangraph, max_c_sources_names, subpangraph.get_nodes_count())
        subsubpangraph = run_poa(subsubpangraph)
        remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)

        subpangraph.pangraph.clear_consensuses()
        subpangraph.pangraph.add_consensus(remapped_best_path)
        # max_c_to_node = subpangraph.get_paths_compatibility_to_consensus(remapped_best_path)
        max_c_to_node = subpangraph.get_paths_compatibility(0)
        remapped_to_orig_best_path = subpangraph.get_consensus_remapped_to_original_nodes(0)
        node_cutoff = find_node_cutoff(max_c_to_node)
        # compatible_sources_ids = get_max_compatible_sources_ids(max_c_to_node, node_cutoff)
        compatible_sources_names = get_max_compatible_sources_ids(current_paths_names, max_c_to_node, node_cutoff)
        # compatible_sources_names = subpangraph.get_sources_names(compatible_sources_ids)

        node = ConsensusNode(sequences_names=list(compatible_sources_names))
        local_consensus_manager.add_node(node, remapped_to_orig_best_path)

        current_paths_names = sorted(list((set(current_paths_names) - set(compatible_sources_names))))
        # current_paths_ids = (set(current_paths_ids)) - set(compatible_sources_ids)
        subpangraph = SubPangraph(s.pangraph, current_paths_names, subpangraph.orig_nodes_count)
    return local_consensus_manager


def get_top_consensus(subpangraph: SubPangraph):
    subpangraph_with_consensus = run_poa(subpangraph)
    return subpangraph_with_consensus.get_consensus_remapped_to_original_nodes(0)


# #todo może jednak consensus manager poza pangraphem powinien być?
# def produce_tree(subpangraph: SubPangraph, consensus_manager: TreeConsensusManager) -> TreeConsensusManager:
#     if node_complete(subpangraph, consensus_manager):
#         return consensus_manager  # ma jeden node i jeden Path w PathManager
#
#     children_cm = get_children_cm(subpangraph)  # consensus manager z węzłami, które są rodzeństwem
#     for child in children_cm.get_nodes():
#         child_subpangraph = SubPangraph(subpangraph.pangraph, child.sequences_ids)
#         child_consensus_manager = produce_tree(child_subpangraph, consensus_manager)
#         remapped_cm = subpangraph.remap_to_original(child_consensus_manager) # każdy consensus zmapować do subpangraph?
#         # child_node_subpangraph = SubPangraph(subpangraph, child.get_sequences_ids())
#         # child_node_subpangraph = subpangraph.keep_paths(child.get_sequences_ids())
#         # child_consensus_manager = produce_tree(child_node_subpangraph, consensus_manager, genomes_info)
#         # remapped_cm = child_node_subpangraph.remap_to_original(child_consensus_manager)
#         consensus_manager.merge(remapped_cm) #uzupełnić listę consensusów (nadać odpowiednie id?), dokleić poddrzewa
#     return consensus_manager


# def get_children_cm(subpangraph: SubPangraph) -> TreeConsensusManager:
#     current_path_ids = subpangraph.get_path_ids()
#     current_path_names = subpangraph.get_sources_names()
#     local_consensus_manager = TreeConsensusManager(max_nodes_count=subpangraph.get_nodes_count())
#     orig_subpangraph = deepcopy(subpangraph)
#     while current_path_names:
#         subpangraph = run_poa(subpangraph)
#         compatibility_to_node_sequences = subpangraph.get_paths_compatibility(0)
#         max_cutoff = find_max_cutoff(compatibility_to_node_sequences)
#         max_compatible_sources_ids = get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff)
#
#         # subsubpangraph = subpangraph.keep_sources_ids(list(max_compatible_sources_ids)) #przemapowanie ale z utratą zupełnie oryginalnego - tu jest niepotrzebny
#         subsubpangraph = SubPangraph(subpangraph.pangraph, list(max_compatible_sources_ids))
#         subsubpangraph = run_poa(subsubpangraph)
#         remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)
#
#         max_compatibility_to_node_sequences = subpangraph.get_paths_compatibility_to_consensus(remapped_best_path)
#         node_cutoff = find_node_cutoff(max_compatibility_to_node_sequences)
#         compatible_sources_ids = get_max_compatible_sources_ids(max_compatibility_to_node_sequences, node_cutoff)
#         compatible_sources_names = subpangraph.get_sources_names(compatible_sources_ids)
#
#         #prace koncowe
#         node = ConsensusNode(sequences_names=list(compatible_sources_names))
#         local_consensus_manager.add_node(node, remapped_best_path)
#
#         current_path_ids = (set(current_path_ids) - set(compatible_sources_ids))
#         # current_path_names = (set(current_path_names) - set(compatible_sources_names))
#         subpangraph = SubPangraph(orig_subpangraph.pangraph, list(current_path_ids), subpangraph.orig_nodes_count) #do kolejnych iteracji tej pętli
#         current_path_names = subpangraph.get_sources_names()
#     return local_consensus_manager


def find_max_cutoff(compatibility_to_node_sequences, cutoff_search_range):
    if not compatibility_to_node_sequences:
        raise ValueError("Empty compatibilities list. Finding max cutoff is not possible.")
    if len(cutoff_search_range) != 2:
        raise ValueError("Cutoff search range must have length 2.")
    if cutoff_search_range[1] < cutoff_search_range[0]:
        raise ValueError("For cutoff search range [x, y] x must be <= y.")

    min_search_pos = round((len(compatibility_to_node_sequences)-1)*cutoff_search_range[0])
    max_search_pos = round((len(compatibility_to_node_sequences)-1)*cutoff_search_range[1])
    sorted_comp = sorted(compatibility_to_node_sequences)
    if min_search_pos == max_search_pos:
        return sorted_comp[min_search_pos]

    search_range = sorted(set(sorted_comp[min_search_pos: max_search_pos+1]))
    if len(search_range) == 1:
        return search_range[0]
    if len(search_range) == 2:
        return search_range[1]

    max_diff = search_range[1] - search_range[0]
    max_cutoff = search_range[0]
    for i in range(1, len(search_range)-1):
        current_diff = search_range[i+1] - search_range[i]
        if current_diff >= max_diff:
            max_diff = current_diff
            max_cutoff = search_range[i+1]

    return max_cutoff


def find_node_cutoff(compatibility_to_node_sequences, multiplier):
    if not compatibility_to_node_sequences:
        raise ValueError("Empty compatibilities list. Finding max cutoff.")
    sorted_comp = sorted(set(compatibility_to_node_sequences))
    if len(compatibility_to_node_sequences) == 1:
        return sorted_comp[0]
    elif len(compatibility_to_node_sequences) == 2:
        return sorted_comp[1]

    mean_distance = (compatibility_to_node_sequences[-1] - compatibility_to_node_sequences[0])/(len(compatibility_to_node_sequences)-1)
    required_gap = mean_distance * multiplier

    distances = np.array([sorted_comp[i + 1] - sorted_comp[i] for i in range(len(sorted_comp)-1)])
    if any(distances >= required_gap):
        a = np.where(distances >= required_gap)[0][0]+1
        return sorted_comp[a]
    else:
        # logger.warning("Cannot find node cutoff for given multiplier. Multiplier == 1 was used instead.")
        return sorted_comp[np.where(distances >= mean_distance)[0][0]+1]


def get_max_compatible_sources_ids(current_paths_names, compatibility_to_node_sequences, max_cutoff):
    npver = np.array(compatibility_to_node_sequences)
    path_names=np.array(current_paths_names)
    return list(path_names[np.where(npver >= max_cutoff)[0]])


def run_poa(subpangraph: SubPangraph) -> SubPangraph:
    pangraph_with_consensus = simple.run(ap.outputdir,
                                         subpangraph.pangraph,
                                         ap.config.hbmin,
                                         ap.genomes_info)
    subpangraph.set_pangraph(pangraph_with_consensus)
    return subpangraph
