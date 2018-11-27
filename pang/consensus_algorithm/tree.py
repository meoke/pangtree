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
from collections import deque

import logging

ap = AlgorithmParams()


def run(outputdir: Path, pangraph: Pangraph, config: TreeConfig, genomes_info: MultialignmentMetadata) -> Pangraph:
    ap.outputdir, ap.config, ap.genomes_info = outputdir, config, genomes_info
    root_consensus_manager = produce_tree(pangraph)
    pangraph.set_consensus_manager(root_consensus_manager)
    return pangraph


def node_ready(node: ConsensusNode):
    min_own_comp = min(node.get_compatibilities_to_own_sources())
    if len(node.sequences_names) in [0, 1] or min_own_comp >= ap.config.stop:
        return True
    return False


def get_root_node(pangraph: Pangraph):
    root_pangraph = SubPangraph(pangraph, pangraph.get_path_names())
    root_node_consensus = get_top_consensus(root_pangraph)
    root_node = ConsensusNode(sequences_names=list(pangraph.get_path_names()))
    root_node.compatibilities_to_all = pangraph.get_paths_compatibility_to_consensus(root_node_consensus)
    root_node.mincomp = min([c for seq, c in root_node.compatibilities_to_all.items() if seq in root_node.sequences_names])
    return root_node, root_node_consensus


def produce_tree(pangraph: Pangraph) -> TreeConsensusManager:
    cm = TreeConsensusManager(max_nodes_count=pangraph.get_nodes_count())
    root_node, root_consensus = get_root_node(pangraph)
    cm.add_node(root_node, root_consensus)
    nodes_to_process = deque([root_node])
    p = deepcopy(pangraph)
    while nodes_to_process:
        subtree_root = nodes_to_process.pop()
        children_nodes_manager = get_children_nodes(p, subtree_root)

        if ap.config.re_consensus:
            children_nodes_manager = reorder_consensuses(p, children_nodes_manager)

        chidren_nodes = children_nodes_manager.get_nodes()
        if len(chidren_nodes) == 1:
            continue

        for child in chidren_nodes:
            consensus = children_nodes_manager.get_consensus(child.consensus_id)
            child.parent_node_id = subtree_root.consensus_id
            child.compatibilities_to_all = pangraph.get_paths_compatibility_to_consensus(consensus)
            child_node_id = cm.add_node(child, consensus)
            subtree_root.children_nodes.append(child_node_id)
            if not node_ready(child):
                nodes_to_process.append(child)

    return cm


def get_children_nodes(orig_pangraph: Pangraph, cn: ConsensusNode) -> TreeConsensusManager:
    current_paths_names = cn.sequences_names
    op = deepcopy(orig_pangraph)
    subpangraph = SubPangraph(op, cn.sequences_names)
    local_consensus_manager = TreeConsensusManager(max_nodes_count=subpangraph.orig_nodes_count)
    logging.info(f"Searching children for id: {cn.consensus_id}, len: {len(cn.sequences_names)}, names: {cn.sequences_names}")
    while current_paths_names:
        subpangraph = run_poa(subpangraph)
        c_to_node = subpangraph.get_paths_compatibility(0)
        max_cutoff = find_max_cutoff(c_to_node, ap.config.cutoff_search_range)
        max_c_sources_names = get_max_compatible_sources_ids(current_paths_names, c_to_node, max_cutoff)

        subsubpangraph = SubPangraph(subpangraph.pangraph, max_c_sources_names, subpangraph.get_nodes_count())
        subsubpangraph = run_poa(subsubpangraph)
        remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)

        subpangraph.pangraph.clear_consensuses()
        subpangraph.pangraph.add_consensus(remapped_best_path)
        max_c_to_node = subpangraph.get_paths_compatibility(0)
        remapped_to_orig_best_path = subpangraph.get_consensus_remapped_to_original_nodes(0)
        if ap.config.anti_granular:
            node_cutoff = find_node_cutoff(max_c_to_node, ap.config.multiplier, local_consensus_manager.get_all_leaves_mincomps())
        else:
            node_cutoff = find_node_cutoff_old(max_c_to_node, ap.config.multiplier)
        compatible_sources_names = get_max_compatible_sources_ids(current_paths_names, max_c_to_node, node_cutoff)

        node = ConsensusNode(sequences_names=list(compatible_sources_names), mincomp=node_cutoff)
        local_consensus_manager.add_node(node, remapped_to_orig_best_path)

        current_paths_names = sorted(list((set(current_paths_names) - set(compatible_sources_names))))
        subpangraph = SubPangraph(op, current_paths_names, subpangraph.orig_nodes_count)

    return local_consensus_manager


def reorder_consensuses(pangraph, tcm: TreeConsensusManager):
    for seqname in tcm.get_sequences_names():
        path = pangraph.get_path(seqname)
        consensus_id_to_comp = {}
        for n in tcm.get_nodes():
            consensus = tcm.get_consensus(n.consensus_id)
            consensus_id_to_comp[n.consensus_id] = pangraph.get_path_compatibility(path, consensus)
            if seqname in n.sequences_names:
                current_consensus_id = n.consensus_id
        best_comp = max(consensus_id_to_comp.values())
        best_comp_id = [comp_id for comp_id, comp_value in consensus_id_to_comp.items()][0]
        if best_comp != consensus_id_to_comp[current_consensus_id]:
            tcm.consensus_tree.nodes[best_comp_id].sequences_names.append(seqname)
            tcm.consensus_tree.nodes[current_consensus_id].sequences_names.remove(seqname)
    for n in tcm.get_nodes():
        if not n.sequences_names:
            tcm.remove_consensus(n.consensus_id)
    return tcm


def get_top_consensus(subpangraph: SubPangraph):
    subpangraph_with_consensus = run_poa(subpangraph)
    return subpangraph_with_consensus.get_consensus_remapped_to_original_nodes(0)


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
    max_cutoff = search_range[1]
    for i in range(1, len(search_range)-1):
        current_diff = search_range[i+1] - search_range[i]
        if current_diff >= max_diff:
            max_diff = current_diff
            max_cutoff = search_range[i+1]

    return max_cutoff


def find_node_cutoff_old(compatibility_to_node_sequences, multiplier):
    if not compatibility_to_node_sequences:
        raise ValueError("Empty compatibilities list. Finding max cutoff.")
    sorted_comp = sorted(set(compatibility_to_node_sequences))
    if len(sorted_comp) == 1:
        return sorted_comp[0]
    elif len(sorted_comp) == 2:
        return sorted_comp[1]

    mean_distance = (sorted_comp[-1] - sorted_comp[0])/(len(sorted_comp)-1)
    required_gap = mean_distance * multiplier

    distances = np.array([sorted_comp[i + 1] - sorted_comp[i] for i in range(len(sorted_comp)-1)])
    if any(distances >= required_gap):
        a = np.where(distances >= required_gap)[0][0]+1
        return sorted_comp[a]
    else:
        logging.warning("Cannot find node cutoff for given multiplier. Multiplier == 1 was used instead.")
        return sorted_comp[np.where(distances >= mean_distance)[0][0]+1]

def find_node_cutoff(compatibility_to_node_sequences, multiplier, mincomps=None):
    if not mincomps:
        return find_node_cutoff_old(compatibility_to_node_sequences, multiplier)
    sorted_comp = sorted(set(compatibility_to_node_sequences))

    if len(sorted_comp) == 1:
        return sorted_comp[0]
    # elif len(sorted_comp) == 2:
    #     return sorted_comp[1]

    anti_granular_guard = min(mincomps) if mincomps else []
    mean_distance = (sorted_comp[-1] - sorted_comp[0])/(len(sorted_comp)-1)
    required_gap = mean_distance * multiplier
    small_comps = [sc for sc in sorted_comp if sc <= anti_granular_guard]

    distances = np.array([small_comps[i + 1] - small_comps[i] for i in range(len(small_comps) - 1)])
    if any(distances >= required_gap):
        a = np.where(distances >= required_gap)[0][0]+1
        return small_comps[a]
    else:
        try:
            left_to_guard = [c for c in sorted_comp if c < anti_granular_guard][-1]
            left_to_guard_diff = anti_granular_guard - left_to_guard
        except:
            left_to_guard_diff = 1
        try:
            right_to_guard = [c for c in sorted_comp if c > anti_granular_guard][0]
            right_to_guard_diff = right_to_guard - anti_granular_guard
        except:
            right_to_guard_diff = 1
        if right_to_guard_diff <= left_to_guard_diff:
            try:
                return right_to_guard
            except:
                pass
        else:
            return left_to_guard

            ###
            #     logging.warning("Cannot find node cutoff for given multiplier. Multiplier == 1 was used instead.")
            #     # return sorted_comp[np.where(distances >= mean_distance)[0][0]+1]


def get_max_compatible_sources_ids(current_paths_names, compatibility_to_node_sequences, max_cutoff):
    npver = np.array(compatibility_to_node_sequences)
    path_names = np.array(current_paths_names)
    return list(path_names[np.where(npver >= max_cutoff)[0]])


def run_poa(subpangraph: SubPangraph) -> SubPangraph:
    pangraph_with_consensus = simple.run(ap.outputdir,
                                         subpangraph.pangraph,
                                         ap.config.hbmin,
                                         ap.genomes_info)
    subpangraph.set_pangraph(pangraph_with_consensus)
    return subpangraph
