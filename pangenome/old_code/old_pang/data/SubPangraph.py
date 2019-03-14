from typing import List
from old_pang.graph.Pangraph import Pangraph
from old_pang.graph.PathManager import PathManager
from old_pang.graph.Node import Node
import numpy as np
import copy


class SubPangraph(object):
    def __init__(self, pangraph: Pangraph, sequences_names_to_keep: List[int], orig_nodes_count: int = None):
        self.pangraph = Pangraph()
        if orig_nodes_count:
            self.orig_nodes_count = orig_nodes_count
        else:
            self.orig_nodes_count = pangraph.get_nodes_count()

        sorted_sequences_names = sorted(sequences_names_to_keep)
        if pangraph.get_path_names() == sorted_sequences_names:
            self.pangraph = copy.deepcopy(pangraph)
            self.nodes_ids_mapping = {i: i for i in range(self.pangraph.get_nodes_count())}
        else:
            self.pangraph._pathmanager = copy.deepcopy(pangraph._pathmanager)
            self.pangraph._pathmanager.keep_paths_names(sorted_sequences_names)
            nodes_ids_to_keep = self.pangraph._pathmanager.get_active_nodes()
            self.pangraph._pathmanager.keep_nodes_ids(nodes_ids_to_keep)
            self.pangraph._nodes, self.nodes_ids_mapping = self.build_nodes(pangraph, nodes_ids_to_keep)
            self.pangraph._consensusmanager = PathManager()

    def build_nodes(self, pangraph: Pangraph, nodes_ids_to_keep: List[int]):
        def is_active_edge(pangraph: Pangraph , node_id: int, in_node: int):
            return True
        new_nodes = [None] * len(nodes_ids_to_keep)
        new_to_old_mapping = {}
        old_to_new_mapping = {}
        for i, node_id in enumerate(nodes_ids_to_keep):
            node = pangraph.get_node(node_id)
            in_nodes = [in_node for in_node in node.in_nodes if in_node in nodes_ids_to_keep]
            in_nodes = [in_node for in_node in in_nodes if is_active_edge(pangraph, node_id, in_node)]
            # todo perf did not work # in_nodes = sorted(set(node.in_nodes) & set(nodes_ids_to_keep))
            next_node_id = node.aligned_to
            while True:
                if next_node_id == node.id or next_node_id is None:
                    aligned_to = None
                    break
                if next_node_id in nodes_ids_to_keep:
                    aligned_to = next_node_id
                    break
                next_node_id = pangraph.get_node(next_node_id).aligned_to
            new_nodes[i] = Node(id=i, base=node.base, in_nodes=in_nodes, aligned_to=aligned_to, column_id=node.column_id, block_id=node.block_id)

            new_to_old_mapping[i] = node_id
            old_to_new_mapping[node_id] = i
        for node in new_nodes:
            if node.aligned_to:
                node.aligned_to = old_to_new_mapping[node.aligned_to]
            node.in_nodes = [old_to_new_mapping[in_node] for in_node in node.in_nodes if in_node in nodes_ids_to_keep]
            # todo perf did not work #node.in_nodes = [old_to_new_mapping[in_node] for in_node in sorted(set(node.in_nodes) & set(nodes_ids_to_keep))]

        return new_nodes, new_to_old_mapping

    def remap_to_original(self, child_consensus_manager):
        raise NotImplemented

    def get_consensus_remapped_to_original_nodes(self, consensus_id):
        original_path = np.zeros(shape=self.orig_nodes_count, dtype=bool)
        consensus_nodes_ids = [self.nodes_ids_mapping[node_id] for node_id in np.where(self.pangraph._consensusmanager.paths[consensus_id])[0]]
        original_path[consensus_nodes_ids] = True
        return original_path

    def get_path_ids(self):
        return self.pangraph.get_path_ids()

    def get_paths_names(self):
        return self.pangraph.get_path_names()

    def get_paths_compatibility(self, consensus_id):
        return self.pangraph.get_paths_compatibility(consensus_id)

    def get_paths_compatibility_to_consensus(self, consensus):
        return self.pangraph.get_paths_compatibility_to_consensus(consensus)

    def set_pangraph(self, pangraph_with_consensus):
        self.pangraph = pangraph_with_consensus

    def keep_sources_ids(self, sources_ids_to_keep):
        return SubPangraph(self.pangraph, sources_ids_to_keep, self.orig_nodes_count)

    def __eq__(self, other):
        return (self.pangraph == other.pangraph and
                self.nodes_ids_mapping == other.nodes_ids_mapping)

    # def get_sources_names(self, specific_sources_ids = None):
    #     if specific_sources_ids is None:
    #         return self.pangraph._pathmanager.get_path_names()
    #     try:
    #         return [self.pangraph.get_source_name(src_id) for src_id in specific_sources_ids]
    #     except NoPath:
    #         return []

    def get_nodes_count(self):
        return self.pangraph.get_nodes_count()



