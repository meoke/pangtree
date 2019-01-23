from typing import List, Dict
from Pangenome import Pangenome


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
    def __init__(self,
                 id: int,
                 genbankID: str,
                 assemblyID: str,
                 mafname: str,
                 name:str,
                 group: str,
                 nodes_ids: List[int]):
        self.id = id
        self.genbankID = genbankID
        self.assemblyID = assemblyID
        self.mafname = mafname
        self.name = name
        self.group = group
        self.nodes_ids = nodes_ids


class JSONConsensus:
    def __init__(self,
                 id: int,
                 name: str,
                 parent: int,
                 children: List[int],
                 comp_to_all_sequences: Dict[str, float],
                 sequences_ids: List[int],
                 nodes_ids: List[int],
                 mincomp: float):
        self.id = id
        self.name = name
        self.parent = parent
        self.children = children
        self.comp_to_all_sequences = comp_to_all_sequences
        self.sequences_ids = sequences_ids
        self.nodes_ids = nodes_ids
        self.mincomp = mincomp


class JSONPangenome:
    def __init__(self, pangenome: Pangenome =None):
        if not pangenome:
            return
        # todo perf # self.nodes = [JSONNode(node.id, decode(node.base)) for node in pangenome.pangraph.get_nodes()]
        self.nodes = []
        self.edges = []
        seqeuences_metadata = [pangenome.genomes_info.genomes_metadata[seqID] for seqID in pangenome.pangraph.get_path_names()]
        self.sequences = []
        pm = pangenome.pangraph._pathmanager
        for i, seq_metadata in enumerate(seqeuences_metadata):
            jsonsequence = JSONSequence(pm.get_path_id(seq_metadata.mafname),
                                       seq_metadata.genbankID,
                                       seq_metadata.assemblyID,
                                       seq_metadata.mafname,
                                       seq_metadata.name,
                                       seq_metadata.group,
                                       []
                          )
            self.sequences.append(jsonsequence)

        cm = pangenome.pangraph._consensusmanager
        cm_tree_nodes = pangenome.pangraph._consensusmanager.consensus_tree.nodes

        self.consensuses = [JSONConsensus(id=node.consensus_id,
                                          name=cm.get_path_name(node.consensus_id),
                                          parent=node.parent_node_id,
                                          children=node.children_nodes,
                                          comp_to_all_sequences={seq_name: float(comp)
                                                                 for seq_name, comp in node.compatibilities_to_all.items()},
                                          sequences_ids=[pm.get_path_id(seq_name)
                                                         for seq_name in node.sequences_names],
                                          # todo perf # nodes_ids = [int(node_id) for node_id in pangenome.pangraph.get_consensus_nodes_ids(cm.get_path_name(node.consensus_id))],
                                          nodes_ids=[],
                                          mincomp=float(node.mincomp)
                                          )
                            for node in cm_tree_nodes]

    # def build_from_dict(self, dictionary):
    #     self.nodes = [JSONNode(node['id'], node['nucleobase']) for node in dictionary['nodes']]
    #     self.edges = dictionary['edges']
    #     self.sequences = [JSONSequence(sequence['id'],
    #                                    sequence['name'],
    #                                    sequence['title'],
    #                                    [int(node_id) for node_id in sequence['nodes_ids']]
    #                                    )
    #                       for sequence in dictionary['sequences']]
    #     self.consensuses = [JSONConsensus(id=n['id'],
    #                                          name=n['name'],
    #                                          parent=n['parent'],
    #                                          children=[c for c in n['children']],
    #                                          comp_to_all_sequences={seqname: comp for seqname, comp in n['comp_to_all_sequences'].items()},
    #                                          sequences_ids=[s for s in n['sequences_ids']],
    #                                          nodes_ids=[c for c in n['nodes_ids']],
    #                                          mincomp=n['mincomp'])
    #                            for n in dictionary['consensuses']]
