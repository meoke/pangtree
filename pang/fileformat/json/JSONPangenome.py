from graph.nucleotides import decode as decode
from typing import List, Dict
import Pangenome


class JSONNode:
    def __init__(self, id: int, nucleobase: str):
        self.id = id
        self.nucleobase = nucleobase


class JSONEdge:
    def __init__(self, source: int, target: int, type: str):
        self.source = source
        self.target = target
        self.type = type


class JSONSequence:
    def __init__(self, id: int, nodes_ids: List[int]):
        self.id = id
        self.nodes_ids = nodes_ids


class JSONConsensus:
    def __init__(self,
                 id: int,
                 name: str,
                 parent: int,
                 children: List[int],
                 comp_to_all_sequences: Dict[str, float],
                 sequences: List[int],
                 nodes_ids: List[int]):
        self.id = id
        self.name = name
        self.parent = parent
        self.children = children
        self.comp_to_all_sequences = comp_to_all_sequences
        self.sequences = sequences
        self.nodes_ids = nodes_ids


class JSONPangenome:
    def __init__(self, pangenome: Pangenome):
        self.nodes = [JSONNode(node.id, decode(node.base)) for node in pangenome.pangraph.get_nodes()]
        self.edges = []
        self.sequences = [JSONSequence(pangenome.genomes_info.get_id(sequence),
                                       [int(node_id) for node_id in pangenome.pangraph.get_sequence_nodes_ids(sequence)])
                          for sequence in pangenome.pangraph.get_path_names()]
        cm = pangenome.pangraph._consensusmanager
        cm_tree_nodes = pangenome.pangraph._consensusmanager.consensus_tree.nodes
        self.consensus_tree = [JSONConsensus(id=node.consensus_id,
                                             name=cm.get_path_name(node.consensus_id),
                                             parent=node.parent_node_id,
                                             children=node.children_nodes,
                                             comp_to_all_sequences={str(seq_id): float(comp)
                               for seq_id, comp in node.compatibilities_to_all.items()} if node.compatibilities_to_all else None,
                                             sequences=node.sequences_names,
                                             nodes_ids=[int(node_id) for node_id in
                                                        pangenome.pangraph.get_consensus_nodes_ids(cm.get_path_name(node.consensus_id))]) for node in cm_tree_nodes]
