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

from collections import deque


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

def produce_tree2(pangraph: Pangraph) -> TreeConsensusManager:
    all_sequences_names = pangraph.get_path_names()  #

    cm = TreeConsensusManager(max_nodes_count=pangraph.get_nodes_count())  # jego będę produkować
    root_node = ConsensusNode(sequences_names=list(all_sequences_names))
    root_pangraph = SubPangraph(pangraph, pangraph.get_path_ids())
    cm.add_node(root_node, get_top_consensus(root_pangraph))

    nodes_to_process = deque([root_node])
    while nodes_to_process:
        subtree_root = nodes_to_process.pop()  # ConsensusTreeNode
        sequences_ids = [pangraph.get_path_id(pathname) for pathname in subtree_root.sequences_names]
        current_node_pangraph = SubPangraph(pangraph, sequences_ids)
        children_nodes = get_children_cm2(current_node_pangraph, subtree_root)  # TreeConsensusManager z węzłami-braćmi i odpowiadającymi im consensusami w wersji dla pangraohu ograniczonego do ścieżek w danym nodzie

        for child in children_nodes.get_nodes():
            consensus_in_subpangraph = children_nodes.get_consensus([child.consensus_id])
            consensus = current_node_pangraph.remap_to_original(consensus_in_subpangraph)

            child.parent_node_id = subtree_root.consensus_id
            child_node_id = cm.add_node(child, consensus)
            subtree_root.children_nodes.append(child_node_id)

            if not node_ready(subtree_root, child_node_id):
                nodes_to_process.append(child)

def get_children_cm2(subpangraph: SubPangraph, node: ConsensusNode) -> TreeConsensusManager:
    current_paths_names = node.sequences_names
    current_paths_ids = [subpangraph.pangraph.get_path_id(path_name) for path_name in current_paths_names]
    local_consensus_manager = TreeConsensusManager(max_nodes_count=subpangraph.get_nodes_count())
    while current_paths_names:
        subpangraph = run_poa(subpangraph)
        c_to_node = subpangraph.get_paths_compatibility(0)
        max_cutoff = find_max_cutoff(c_to_node)
        max_c_sources_ids = get_max_compatible_sources_ids(c_to_node, max_cutoff)

        subsubpangraph = SubPangraph(subpangraph.pangraph, max_c_sources_ids, subpangraph.get_nodes_count())
        subsubpangraph = run_poa(subsubpangraph)
        remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)

        subpangraph.pangraph.clear_consensuses()
        subpangraph.pangraph.add_consensus(remapped_best_path)
        # max_c_to_node = subpangraph.get_paths_compatibility_to_consensus(remapped_best_path)
        max_c_to_node = subpangraph.get_paths_compatibility(0)
        remapped_to_orig_best_path = subpangraph.get_consensus_remapped_to_original_nodes(0)
        node_cutoff = find_node_cutoff(max_c_to_node)
        compatible_sources_ids = get_max_compatible_sources_ids(max_c_to_node, node_cutoff)
        compatible_sources_names = subpangraph.get_sources_names(compatible_sources_ids)

        node = ConsensusNode(sequences_names=list(compatible_sources_names))
        local_consensus_manager.add_node(node, remapped_to_orig_best_path)

        current_paths_names = (set(current_paths_names) - set(compatible_sources_names))
        subpangraph = SubPangraph(subpangraph.pangraph, current_paths_ids, subpangraph.orig_nodes_count)
    return local_consensus_manager

def get_top_consensus(subpangraph: SubPangraph):
    subpangraph_with_consensus = run_poa(subpangraph)
    return subpangraph_with_consensus.get_consensus_remapped_to_original_nodes(0)


#todo może jednak consensus manager poza pangraphem powinien być?
def produce_tree(subpangraph: SubPangraph, consensus_manager: TreeConsensusManager) -> TreeConsensusManager:
    if node_complete(subpangraph, consensus_manager):
        return consensus_manager  # ma jeden node i jeden Path w PathManager

    children_cm = get_children_cm(subpangraph)  # consensus manager z węzłami, które są rodzeństwem
    for child in children_cm.get_nodes():
        child_subpangraph = SubPangraph(subpangraph.pangraph, child.sequences_ids)
        child_consensus_manager = produce_tree(child_subpangraph, consensus_manager)
        remapped_cm = subpangraph.remap_to_original(child_consensus_manager) # każdy consensus zmapować do subpangraph?
        # child_node_subpangraph = SubPangraph(subpangraph, child.get_sequences_ids())
        # child_node_subpangraph = subpangraph.keep_paths(child.get_sequences_ids())
        # child_consensus_manager = produce_tree(child_node_subpangraph, consensus_manager, genomes_info)
        # remapped_cm = child_node_subpangraph.remap_to_original(child_consensus_manager)
        consensus_manager.merge(remapped_cm) #uzupełnić listę consensusów (nadać odpowiednie id?), dokleić poddrzewa
    return consensus_manager


def get_children_cm(subpangraph: SubPangraph) -> TreeConsensusManager:
    current_path_ids = subpangraph.get_path_ids()
    current_path_names = subpangraph.get_sources_names()
    local_consensus_manager = TreeConsensusManager(max_nodes_count=subpangraph.get_nodes_count())
    orig_subpangraph = deepcopy(subpangraph)
    while current_path_names:
        subpangraph = run_poa(subpangraph)
        compatibility_to_node_sequences = subpangraph.get_paths_compatibility(0)
        max_cutoff = find_max_cutoff(compatibility_to_node_sequences)
        max_compatible_sources_ids = get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff)

        # subsubpangraph = subpangraph.keep_sources_ids(list(max_compatible_sources_ids)) #przemapowanie ale z utratą zupełnie oryginalnego - tu jest niepotrzebny
        subsubpangraph = SubPangraph(subpangraph.pangraph, list(max_compatible_sources_ids))
        subsubpangraph = run_poa(subsubpangraph)
        remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)

        max_compatibility_to_node_sequences = subpangraph.get_paths_compatibility_to_consensus(remapped_best_path)
        node_cutoff = find_node_cutoff(max_compatibility_to_node_sequences)
        compatible_sources_ids = get_max_compatible_sources_ids(max_compatibility_to_node_sequences, node_cutoff)
        compatible_sources_names = subpangraph.get_sources_names(compatible_sources_ids)

        #prace koncowe
        node = ConsensusNode(sequences_names=list(compatible_sources_names))
        local_consensus_manager.add_node(node, remapped_best_path)

        current_path_ids = (set(current_path_ids) - set(compatible_sources_ids))
        # current_path_names = (set(current_path_names) - set(compatible_sources_names))
        subpangraph = SubPangraph(orig_subpangraph.pangraph, list(current_path_ids), subpangraph.orig_nodes_count) #do kolejnych iteracji tej pętli
        current_path_names = subpangraph.get_sources_names()
    return local_consensus_manager


def find_max_cutoff(compatibility_to_node_sequences):
    cutoff_search_range = ap.config.cutoff_search_range
    #określić które compatiilities przeglądać -> indeksy
    min_search_idx = round((len(compatibility_to_node_sequences))*cutoff_search_range[0])
    max_search_idx = round((len(compatibility_to_node_sequences))*cutoff_search_range[1])
    #znaleźć wśród nich największą różnicę
    compatibilities_to_be_searched = sorted(compatibility_to_node_sequences)[min_search_idx: max_search_idx]
    if len(compatibilities_to_be_searched) == 1 or len(compatibilities_to_be_searched) == 2:
        return compatibilities_to_be_searched[0]
    elif len(compatibilities_to_be_searched) == 0:
        return sorted(compatibility_to_node_sequences)[min_search_idx]
    max_diff = 0
    differences = []
    for i, c in enumerate(compatibilities_to_be_searched):
        if i < len(compatibilities_to_be_searched)-1:
            differences.append(compatibilities_to_be_searched[i+1] - compatibilities_to_be_searched[i])
    max_difference = max(differences)
    cutoff = compatibilities_to_be_searched[differences.index(max_difference)+1]
    return cutoff


def find_node_cutoff(compatibility_to_node_sequences):
    #todo naprawić
    return 0.9


def get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff):
    return np.where(np.array(compatibility_to_node_sequences) >= max_cutoff)[0]


def node_complete(subpangraph, consensus_manager):
    node = consensus_manager.get_root_node()
    if len(node.sequences_ids) == 1:
        return True
    try:
        consensus = consensus_manager.get_consensus(node.consensus_id)
        node_compatibility = min(subpangraph.get_paths_compatibility_to_consensus(consensus))  # czy nodes ids się zgadzają?
        if node_compatibility >= ap.config.stop:
            return True
    except NoConsensus:
        return False



def run_poa(subpangraph: SubPangraph) -> SubPangraph:
    pangraph_with_consensus = simple.run(ap.outputdir,
                                         subpangraph.pangraph,
                                         ap.config.hbmin,
                                         ap.genomes_info)
    subpangraph.set_pangraph(pangraph_with_consensus)
    return subpangraph
