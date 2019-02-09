from io import StringIO
from pathlib import Path
from typing import List, Dict
from fileformats.maf.DAGMaf import DAGMaf
from graph.PangraphBuilder.PangraphBuilderBase import PangraphBuilderBase
from graph.PangraphBuilder.PangraphBuilderFromDAG import PangraphBuilderFromDAG
from graph.PangraphBuilder.PangraphBuilderFromMAF import PangraphBuilderFromMAF
from .Node import Node
from .PathManager import PathManager
import numpy as np
from .FastaSource import EntrezFastaSource, FastaFileSystemSource


class Pangraph:
    def __init__(self):
        self._nodes = []
        self._pathmanager = PathManager()
        self._consensusmanager = PathManager()

    # def __init__(self, max_nodes_count: int=0, start_node_id: int=0, paths_names: List[str]=None):
    #     self._nodes = [None] * max_nodes_count
    #     self._pathmanager = PathManager(start_node_id, max_nodes_count, paths_names)
    #     self._consensusmanager = PathManager()

    def __eq__(self, other):
        return (self._nodes == other._nodes and
                self._pathmanager == other._pathmanager and
                self._consensusmanager == other._consensusmanager)

    def build_from_maf_firstly_converted_to_dag(self, mafcontent: StringIO, fasta_source, genomes_info):
        builder: PangraphBuilderBase = PangraphBuilderFromDAG(genomes_info, fasta_source)
        self._build(builder, mafcontent)

    def build_from_maf(self, mafcontent, genomes_info):
        builder: PangraphBuilderBase = PangraphBuilderFromMAF(genomes_info)
        self._build(builder, mafcontent)

    def _build(self, builder, build_input):
        builder.build(build_input, self)

    # def update(self, pangraph, start_node):
    #     self.update_nodes(pangraph._nodes)
    #     self._pathmanager.update(pangraph._pathmanager, start=start_node)

    def update_nodes(self, new_nodes: List[Node]):
        #todo something to control new_nodes correctness
        if not new_nodes:
            raise Exception("empty new nodes")
        if len(self._nodes) <= new_nodes[-1].id:
            self._nodes = new_nodes
            return
        self._nodes[new_nodes[0].id: new_nodes[-1].id] = new_nodes

    def get_nodes_count(self):
        return len(self._nodes)

    def get_nodes(self):
        return self._nodes

    # def trim_nodes(self, nodes_count: int):
    #     del self._nodes[nodes_count:]
    #     self._pathmanager.trim(nodes_count)

    def set_paths(self, max_nodes_count: int, paths_to_node_ids: Dict[str, List[int]] = None):
        # todo something to control paths correctness
        self._pathmanager.init_from_dict(max_nodes_count, paths_to_node_ids)

    def add_path_to_node(self, path_name, node_id):
        self._pathmanager.mark_and_add(path_name, node_id)

    def get_in_nodes(self, node_id):
        return self._pathmanager.get_in_nodes(node_id)

    # def add_node(self, node: Node, node_id: str):
    #     self._nodes[node_id] = node

    def fill_in_nodes(self):
        for node in self._nodes:
            node.in_nodes = self.get_in_nodes(node.id)

    def get_paths_count(self):
        return self._pathmanager.get_paths_count()

    def get_path_names(self):
        return self._pathmanager.get_path_names()

    def get_path_ids(self):
        return self._pathmanager.get_path_ids()

    def get_path_id(self, pathname):
        return self._pathmanager.get_path_id(pathname)

    def get_path_nodes_count(self, pathname):
        return self._pathmanager.get_path_nodes_count(pathname)

    def get_start_node_id(self, source):
        return self._pathmanager.get_start_node_id(source)

    def get_sources_weights_dict(self):
        return self._pathmanager.get_sources_weights_dict()

    def get_source_consensus_id(self, source):
        return -1

    def get_sources_ids(self, node_id: int) -> List[int]:
        return self._pathmanager.get_sources_ids(node_id)

    def set_consensuses(self, max_nodes_count, paths_to_node_ids):
        self._consensusmanager.init_from_dict(max_nodes_count, paths_to_node_ids)

    def set_cm(self, cm):
        self._consensusmanager = cm

    def get_top_consensus(self):
        return self._consensusmanager.get_top_consensus()

    def get_node(self, node_id):
        return self._nodes[node_id]

    def remove_empty_paths(self):
        self._pathmanager.remove_empty_paths()

    def get_paths(self):
        return self._pathmanager.get_paths()

    def get_path(self, pathname):
        return self._pathmanager.get_path(pathname)

    def get_path_compatibility(self, path, consensus):
        common_nodes_count = np.sum(path & consensus)
        source_nodes_count = np.sum(path)
        return round(common_nodes_count / source_nodes_count, 3)
        # return common_nodes_count / source_nodes_count

    def get_paths_compatibility(self, consensus_id):
        consensus = self._consensusmanager.paths[consensus_id]
        return [self.get_path_compatibility(path, consensus) for path in self._pathmanager.paths]

    def get_paths_compatibility_to_consensus(self, consensus):
        return {self._pathmanager.get_path_name(path_id): self.get_path_compatibility(path, consensus)
                for path_id, path in enumerate(self._pathmanager.paths)}

    def add_consensus(self, consensus):
        self._consensusmanager.add_path("CONSENSUS", consensus)

    def get_source_name(self, src_id):
        return self._pathmanager.get_path_name(src_id)

    def clear_consensuses(self):
        self._consensusmanager.clear_paths()

    def set_consensus_manager(self, consensus_manager):
        self._consensusmanager = consensus_manager

    def get_sequence_nodes_ids(self, sequence):
        return self._pathmanager.get_nodes_ids(sequence)

    def get_consensus_nodes_ids(self, sequence):
        return self._consensusmanager.get_nodes_ids(sequence)
