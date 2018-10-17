from typing import List
from graph.Pangraph import Pangraph
from graph.PathManager import PathManager
from graph.Node import Node


class SubPangraph(object):
    def __init__(self, pangraph: Pangraph, sequences_ids_to_keep: List[int]):
        self.pangraph = Pangraph()

        if pangraph.get_path_ids() == sequences_ids_to_keep:
            self.pangraph = pangraph
            self.nodes_ids_mapping = {i: i for i in range(self.pangraph.get_nodes_count())}
        else:
            self.pangraph._pathmanager = pangraph._pathmanager
            self.pangraph._pathmanager.keep_paths_ids(sequences_ids_to_keep)
            nodes_ids_to_keep = self.pangraph._pathmanager.get_active_nodes()
            self.pangraph._pathmanager.keep_nodes_ids(nodes_ids_to_keep)
            self.pangraph._nodes, self.nodes_ids_mapping = self.build_nodes(pangraph, nodes_ids_to_keep)
            self.pangraph._consensusmanager = PathManager()

    # def set_pangraph(self, pangraph):
    #     self.pangraph = pangraph

    # def get_mapping(self):
    #     return self.nodes_ids_mapping

    # def keep_paths(self, sequences_ids_to_keep: List[int]):
    #     raise NotImplemented

    def build_nodes(self, pangraph: Pangraph, nodes_ids_to_keep: List[int]):
        new_nodes = [None] * len(nodes_ids_to_keep)
        new_to_old_mapping = {}
        old_to_new_mapping = {}
        for i, node_id in enumerate(nodes_ids_to_keep):
            node = pangraph.get_node(node_id)
            in_nodes = [in_node for in_node in node.in_nodes if in_node in nodes_ids_to_keep]
            next_node_id = node.aligned_to
            while True:
                if next_node_id == node.id or next_node_id is None:
                    aligned_to = None
                    break
                if next_node_id in nodes_ids_to_keep:
                    aligned_to = next_node_id
                    break
                next_node_id = pangraph.get_node(next_node_id).aligned_to
            new_nodes[i] = Node(id=i, base=node.base, in_nodes=in_nodes, aligned_to=aligned_to)

            new_to_old_mapping[i] = node_id
            old_to_new_mapping[node_id] = i
        for node in new_nodes:
            if node.aligned_to:
                node.aligned_to = old_to_new_mapping[node.aligned_to]
            node.in_nodes = [old_to_new_mapping[in_node] for in_node in node.in_nodes if in_node in nodes_ids_to_keep]

        return new_nodes, new_to_old_mapping

    def remap_to_original(self, child_consensus_manager):
        raise NotImplemented

    def get_consensus_remapped_to_original_nodes(self, consensus_id):
        return [self.nodes_ids_mapping[node_id] for node_id in self.pangraph._consensusmanager.paths[consensus_id]]

    def get_path_ids(self):
        return self.pangraph.get_path_ids()

    def get_paths_compatibility(self, consensus_id):
        return self.pangraph.get_paths_compatibility(consensus_id)

    def get_paths_compatibility_to_consensus(self, consensus):
        return self.pangraph.get_paths_compatibility_to_consensus(consensus)

    def set_pangraph(self, pangraph_with_consensus):
        self.pangraph = pangraph_with_consensus

    def keep_sources_ids(self, sources_ids_to_keep):
        return SubPangraph(self.pangraph, sources_ids_to_keep)

    def __eq__(self, other):
        return (self.pangraph == other.pangraph and
                self.nodes_ids_mapping == other.nodes_ids_mapping)