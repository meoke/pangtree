from io import StringIO
from typing import List, Dict
import numpy as np

from graph2.PangraphBuilder.PangraphBuilderBase import PangraphBuilderBase
from graph2.PangraphBuilder.PangraphBuilderFromDAG import PangraphBuilderFromDAG
from graph2.PangraphBuilder.PangraphBuilderFromMAF import PangraphBuilderFromMAF
from graph2.FastaSource import FastaSource
from metadata.MultialignmentMetadata import MultialignmentMetadata
from .Node import Node
from .custom_types import NodeID, SequenceID


class Pangraph:
    def __init__(self):
        self.nodes: List[Node] = []
        self.paths: Dict[SequenceID, List[NodeID]] = {}

    def __eq__(self, other):
        return (self.nodes == other.nodes and
                self.paths == other.paths)

    def build_from_maf_firstly_converted_to_dag(self, mafcontent: StringIO, fasta_source: FastaSource, genomes_info: MultialignmentMetadata):
        builder: PangraphBuilderBase = PangraphBuilderFromDAG(genomes_info, fasta_source)
        self._build(builder, mafcontent)

    def build_from_maf(self, mafcontent: StringIO, genomes_info: MultialignmentMetadata):
        builder: PangraphBuilderBase = PangraphBuilderFromMAF(genomes_info)
        self._build(builder, mafcontent)

    def _build(self, builder: PangraphBuilderBase, build_input):
        builder.build(build_input, self)

    def get_compatibilities(self, sequences_ids: List[SequenceID], consensus: List[NodeID]):
        compatibilities = dict()
        for seq_id in sequences_ids:
            try:
                paths = self.paths[seq_id]
            except KeyError:
                raise Exception("No sequence with given ID in pangraph.")
            if len(paths) == 1:
                path = paths[0]
            else:
                path = [node_id for path in paths for node_id in path]
            compatibilities[seq_id] = len(set(path).intersection(set(consensus)))/len(path)
        return compatibilities

    def get_sequence_nodes_count(self, seq_id):
        if seq_id not in self.paths:
            raise Exception("No sequence with given ID in pangraph.")
        return sum([len(path) for path in self.paths[seq_id]])

    def get_sequences_weights(self, sequences_ids):
        if not sequences_ids:
            return dict()

        a = np.zeros(len(self.nodes), dtype=np.int)
        unweighted_sources_weights = {}
        for seq_id in sequences_ids:
            for path in self.paths[seq_id]:
                a[path] += 1

        for seq_id in sequences_ids:
            sequence_nodes_ids = [node_id for path in self.paths[seq_id] for node_id in path]
            unweighted_sources_weights[seq_id] = np.mean(a[sequence_nodes_ids])

        max_weight = max(unweighted_sources_weights.values())
        min_weight = min(unweighted_sources_weights.values())
        diff_weight = max_weight - min_weight
        if diff_weight == 0:
            normalized_sources_weights_dict = {path_key: 100 for path_key in unweighted_sources_weights.keys()}
        else:
            normalized_sources_weights_dict = {path: int((weight - min_weight)/diff_weight*100) for path, weight in unweighted_sources_weights.items()}
        return normalized_sources_weights_dict

    # def update_nodes(self, new_nodes: List[Node]):
    #     #todo something to control new_nodes correctness
    #     if not new_nodes:
    #         raise Exception("empty new nodes")
    #     if len(self._nodes) <= new_nodes[-1].id:
    #         self._nodes = new_nodes
    #         return
    #     self._nodes[new_nodes[0].id: new_nodes[-1].id] = new_nodes

    # def get_nodes_count(self):
    #     return len(self._nodes)

    # def set_paths(self, max_nodes_count: int, paths_to_node_ids: Dict[str, List[int]] = None):
    #     # todo something to control paths correctness
    #     self._pathmanager.init_from_dict(max_nodes_count, paths_to_node_ids)
    #
    # def add_path_to_node(self, path_name, node_id):
    #     self._pathmanager.mark_and_add(path_name, node_id)

    # def get_in_nodes(self, node_id):
    #     return self._pathmanager.get_in_nodes(node_id)

    # def add_node(self, node: Node, node_id: str):
    #     self._nodes[node_id] = node

    # def fill_in_nodes(self):
    #     for node in self._nodes:
    #         node.in_nodes = self.get_in_nodes(node.id)

    #
    # def get_path_ids(self):
    #     return self._pathmanager.get_path_ids()
    #
    # def get_path_id(self, pathname):
    #     return self._pathmanager.get_path_id(pathname)
    #
    # def get_path_nodes_count(self, pathname):
    #     return self._pathmanager.get_path_nodes_count(pathname)

    # def get_start_node_id(self, source):
    #     return self._pathmanager.get_start_node_id(source)

    # def get_sources_weights_dict(self):
    #     return self._pathmanager.get_sources_weights_dict()

    # def get_source_consensus_id(self, source):
    #     return -1

    # def get_sources_ids(self, node_id: int) -> List[int]:
    #     return self._pathmanager.get_sources_ids(node_id)

    # def set_consensuses(self, max_nodes_count, paths_to_node_ids):
    #     self._consensusmanager.init_from_dict(max_nodes_count, paths_to_node_ids)
    #
    # def get_top_consensus(self):
    #     return self._consensusmanager.get_top_consensus()
    #
    # def get_node(self, node_id):
    #     return self._nodes[node_id]
    #
    # def remove_empty_paths(self):
    #     self._pathmanager.remove_empty_paths()
    #
    # def remove_empty_nodes(self):
    #     nodes_to_remove = []
    #     for i, node in enumerate(self._nodes):
    #         if node is None:
    #             nodes_to_remove.append(i)
    #     for i in sorted(nodes_to_remove, reverse=True):
    #         del self._nodes[i]
    #     if nodes_to_remove:
    #         self._pathmanager.remove_nodes_greater_then(min(nodes_to_remove))
    #
    # def get_paths(self):
    #     return self._pathmanager.get_paths()
    #
    # def get_path(self, pathname):
    #     return self._pathmanager.get_path(pathname)
    #
    # def get_path_compatibility(self, path, consensus):
    #     common_nodes_count = np.sum(path & consensus)
    #     source_nodes_count = np.sum(path)
    #     # return round(common_nodes_count / source_nodes_count, 3)
    #     return common_nodes_count / source_nodes_count
    #
    # def get_paths_compatibility(self, consensus_id):
    #     consensus = self._consensusmanager.paths[consensus_id]
    #     return [self.get_path_compatibility(path, consensus) for path in self._pathmanager.paths]
    #
    # def get_paths_compatibility_to_consensus(self, consensus):
    #     return {self._pathmanager.get_path_name(path_id): self.get_path_compatibility(path, consensus)
    #             for path_id, path in enumerate(self._pathmanager.paths)}
    #
    # def add_consensus(self, consensus):
    #     self._consensusmanager.add_path("CONSENSUS", consensus)
    #
    # def get_source_name(self, src_id):
    #     return self._pathmanager.get_path_name(src_id)
    #
    # def clear_consensuses(self):
    #     self._consensusmanager.clear_paths()
    #
    # def set_consensus_manager(self, consensus_manager):
    #     self._consensusmanager = consensus_manager
    #
    # # def get_sequence_nodes_ids(self, sequence):
    # #     return self._pathmanager.get_nodes_ids(sequence)
    #
    # # def get_consensus_nodes_ids(self, sequence):
    # #     return self._consensusmanager.get_nodes_ids(sequence)
