from graph.Pangraph import Pangraph
from pathlib import Path
from metadata.MultialignmentMetadata import MultialignmentMetadata
from . import simple
import numpy as np
from consensus_algorithm.TreeConfig import TreeConfig
from .AlgorithmParams import AlgorithmParams
from consensus_data.SubPangraph import SubPangraph
from consensus_data.TreeConsensusManager import TreeConsensusManager
from consensus_data.ConsensusNode import ConsensusNode
from consensus_data.Errors import NoConsensus


algorithm_params = AlgorithmParams()


def run(outputdir: Path, pangraph: Pangraph, config: TreeConfig, genomes_info: MultialignmentMetadata) -> Pangraph:
    algorithm_params.outputdir = outputdir
    algorithm_params.config = config
    algorithm_params.genomes_info = genomes_info

    all_sequences_ids = pangraph.get_path_ids()
    subpangraph = SubPangraph(pangraph, all_sequences_ids)
    cm = TreeConsensusManager()
    root_consensusManager = produce_tree(subpangraph, cm) # todo z takimi parametrami, żeby był root z jednym consensusem

    pangraph.set_consensus_manager(root_consensusManager)
    return pangraph

#todo może jednak consensus manager poza pangraphem powinien być?
def produce_tree(subpangraph: SubPangraph, consensus_manager: TreeConsensusManager) -> TreeConsensusManager:
    if node_complete(subpangraph, consensus_manager):
        return consensus_manager
        # return #todo co zwrocic? #consensus manager w którym jest jeden ConsensusNode w ConsensusTree i jeden Path w jego PathManager

    children_cm = get_children_cm(subpangraph)  # consensus manager z węzłami, które są rodzeństwem
    for child in children_cm.get_nodes():
        child_subpangraph = SubPangraph(subpangraph.pangraph, child.get_sequences_ids())  # czy tak mogłoby zostać?
        child_consensus_manager = produce_tree(child_subpangraph, consensus_manager)
        remapped_cm = subpangraph.remap_to_original(child_consensus_manager)
        # child_node_subpangraph = SubPangraph(subpangraph, child.get_sequences_ids())
        # child_node_subpangraph = subpangraph.keep_paths(child.get_sequences_ids())
        # child_consensus_manager = produce_tree(child_node_subpangraph, consensus_manager, genomes_info)
        # remapped_cm = child_node_subpangraph.remap_to_original(child_consensus_manager)
        consensus_manager.merge(remapped_cm) #uzupełnić listę consensusów (nadać odpowiednie id?), dokleić poddrzewa
    return consensus_manager


def get_children_cm(subpangraph) -> TreeConsensusManager:
    current_path_ids = subpangraph.get_path_ids()
    local_consensus_manager = TreeConsensusManager()
    while current_path_ids:
        subpangraph = run_poa(subpangraph)
        compatibility_to_node_sequences = subpangraph.get_paths_compatibility(0)
        max_cutoff = find_max_cutoff(compatibility_to_node_sequences)
        max_compatible_sources_ids = get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff)

        #operacje na subsubpangraphie
        #todo wywalić poza tą funkcję
        subsubpangraph = subpangraph.keep_sources_ids(max_compatible_sources_ids) #przemapowanie ale z utratą zupełnie oryginalnego - tu jest niepotrzebny
        subsubpangraph = run_poa(subsubpangraph)
        remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)

        #cutoff dla wszystkich!!!
        compatibility_to_node_sequences = subpangraph.get_paths_compatibility_to_consensus(remapped_best_path)
        node_cutoff = find_node_cutoff(compatibility_to_node_sequences)
        compatible_sources_ids = get_max_compatible_sources_ids(compatibility_to_node_sequences, node_cutoff)

        #prace koncowe
        subpangraph = subpangraph.keep_sources_ids(compatible_sources_ids) #do kolejnych iteracji tej pętli
        node = ConsensusNode(sequences_ids=compatible_sources_ids)
        local_consensus_manager.add_node(node)
        current_path_ids = (set(current_path_ids) - set(compatible_sources_ids))
    return local_consensus_manager


def find_max_cutoff(compatibility_to_node_sequences):
    cutoff_search_range = algorithm_params.config.cutoff_search_range
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
    pass


def get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff):
    return np.where(compatibility_to_node_sequences >= max_cutoff)[0]


def node_complete(subpangraph, consensus_manager):
    node = consensus_manager.get_root_node()
    if len(node.sequences_ids) == 1:
        return True
    try:
        consensus = consensus_manager.get_consensus(node.consensus_id)
        node_compatibility = subpangraph.get_paths_compatibility_to_consensus(consensus)  # czy nodes ids się zgadzają?
        if node_compatibility >= algorithm_params.config.stop:
            return True
    except NoConsensus:
        return False



def run_poa(subpangraph: SubPangraph) -> SubPangraph:
    pangraph_with_consensus = simple.run(algorithm_params.outputdir,
                                         subpangraph.pangraph,
                                         algorithm_params.config.hbmin,
                                         algorithm_params.genomes_info)
    subpangraph.set_pangraph(pangraph_with_consensus)
    return subpangraph
