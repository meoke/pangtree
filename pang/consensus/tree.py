from graph.Pangraph import Pangraph
import po.writer as powriter
import po.reader as poreader
import userio.pathtools as pathtools
from pathlib import Path
import consensus.poa as poa
from metadata.MultialignmentMetadata import MultialignmentMetadata
from graph.PathManager import PathManager
from graph.Node import Node
from . import simple
from typing import List
import numpy as np
from consensus.TreeConfig import TreeConfig


class ConsensusNode(object):
    def __init__(self):
        self.children_nodes = []
        self.sequences_ids = []
        self.consensus_id = -1


class ConsensusTree(object):
    pass


class ConsensusManager(PathManager):
    def __init__(self):
        self.consensus_tree = ConsensusTree()


class SubPangraph(Pangraph):
    def __init__(self, pangraph: Pangraph, sequences_ids: List[int]):
        Pangraph.__init__(self)
        self._pathmanager = pangraph._pathmanager
        self._pathmanager.keep_paths_ids(sequences_ids)
        nodes_ids_to_keep = self._pathmanager.get_active_nodes()
        self._nodes, self.nodes_ids_mapping = self.build_nodes(pangraph, nodes_ids_to_keep)
        self._consensusmanager = pangraph._consensusmanager

    def build_nodes(self, pangraph: Pangraph, nodes_ids_to_keep: List[int]):
        new_nodes = [None] * len(nodes_ids_to_keep)
        mapping = {}
        for i, node_id in enumerate(nodes_ids_to_keep):
            node = pangraph.get_node(node_id)
            in_nodes = [in_node for in_node in node.in_nodes if in_node in nodes_ids_to_keep]
            while True:
                next_node_id = node.aligned_to
                if next_node_id is None or next_node_id in nodes_ids_to_keep:
                    aligned_to = next_node_id
                    break
                else:
                    next_node_id = pangraph.get_node(next_node_id).aligned_to
            new_nodes[i] = Node(id=i, base=node.base, in_nodes=in_nodes, aligned_to=aligned_to)
            mapping[i] = node_id
        return new_nodes, mapping

    def remap_to_original(self, child_consensus_manager):
        pass

    def get_consensus_remapped_to_original_nodes(self, consensus_id):
        return [self.nodes_ids_mapping[node_id] for node_id in self._consensusmanager.paths[consensus_id]]


def run(outputdir: Path, pangraph: Pangraph, config: TreeConfig, genomes_info: MultialignmentMetadata) -> Pangraph:
    all_sequences_ids = pangraph.get_path_ids()
    subpangraph = SubPangraph(pangraph, all_sequences_ids)  # ciężar przygotowania przemapowań
    cm = ConsensusManager()
    root_consensusManager = produce_tree(outputdir, subpangraph, cm, config, genomes_info) # z takimi parametrami, żeby był root z jednym consensusem

    pangraph.set_consensus_manager(root_consensusManager)
    return pangraph


def produce_tree(outputdir, subpangraph: SubPangraph, consensus_manager: ConsensusManager, config: TreeConfig, genomes_info) -> ConsensusManager:
    if node_complete(subpangraph, consensus_manager):
        return #todo co zwrocic?

    children = get_children(outputdir, subpangraph, config, genomes_info)
    for child in children:
        child_node_subpangraph = SubPangraph(subpangraph, child.get_sequences_ids())
        child_consensus_manager = produce_tree(child_node_subpangraph, consensus_manager, genomes_info)
        remapped_cm = child_node_subpangraph.remap_to_original(child_consensus_manager)
        consensus_manager.merge(remapped_cm)
    return consensus_manager





def find_max_cutoff(compatibility_to_node_sequences, cutoff_search_range):
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


def find_node_cutoff(compatibility_to_node_sequences, multiplier, smallest_comp_up_to_now):
    pass


def get_children(outputdir, subpangraph, config: TreeConfig, genomes_info):
    current_path_ids = subpangraph.get_path_ids()
    children_nodes = []
    while current_path_ids:
        subpangraph = run_poa(outputdir, subpangraph, config.hbmin, genomes_info)
        compatibility_to_node_sequences = subpangraph.get_paths_compatibility(0)
        max_cutoff = find_max_cutoff(compatibility_to_node_sequences, config.cutoff_search_range)
        max_compatible_sources_ids = get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff)

        #operacje na subsubpangraphie
        #todo wywalić poza tą funkcję
        #todo tu jest problem
        subsubpangraph = SubPangraph(subpangraph, max_compatible_sources_ids)
        subsubpangraph = run_poa(outputdir, subsubpangraph, config.hbmin, genomes_info)
        remapped_best_path = subsubpangraph.get_consensus_remapped_to_original_nodes(0)

        #cutoff dla wszystkich!!!
        compatibility_to_node_sequences = subpangraph.get_paths_compatibility(remapped_best_path)
        node_cutoff = find_node_cutoff(compatibility_to_node_sequences, config.multiplier, config.smallest_comp_up_to_now)
        compatible_sources_ids = get_max_compatible_sources_ids(compatibility_to_node_sequences, node_cutoff)

        #prace koncowe
        subpangraph = SubPangraph(subpangraph, compatible_sources_ids)
        tree_node = ConsensusNode(compatible_sources_ids, remapped_best_path)
        children_nodes.append(tree_node)
        current_path_ids = (set(current_path_ids) - set(compatible_sources_ids))



def get_max_compatible_sources_ids(compatibility_to_node_sequences, max_cutoff):
    return np.where(compatibility_to_node_sequences >= max_cutoff)[0]


def node_complete(subpangraph, consensus_manager):
    pass


def run_poa(outputdir, pangraph, hbmin: float, genomes_info):
    pangraph_with_consensus = simple.run(outputdir, pangraph, hbmin, genomes_info)
    return pangraph_with_consensus
