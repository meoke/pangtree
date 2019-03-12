import itertools
from typing import List, Dict
from Pangenome import Pangenome
from arguments.PangenomeParameters import PangenomeParameters, FastaComplementationOption


class JSONProgramParameters:
    def __init__(self, params: PangenomeParameters):
        self.multialignment_file_path: str = str(params.multialignment_file_path)
        self.metadata_file_path: str = str(params.metadata_file_path)
        self.output_path: str = str(params.output_path)
        self.generate_fasta: bool = params.generate_fasta
        self.consensus_type: str = str(params.consensus_type)
        self.max_cutoff_strategy = params.max_cutoff_option
        self.node_cutoff_strategy = params.node_cutoff_option
        self.hbmin: float = params.hbmin
        self.r: float = params.search_range
        self.multiplier: float = params.multiplier
        self.stop: float = params.stop
        self.re_consensus: bool = params.re_consensus
        self.not_dag: bool = params.not_dag
        self.fasta_complementation_option: FastaComplementationOption = str(params.fasta_complementation_option)
        self.local_fasta_dirpath: str = str(params.local_fasta_dirpath)


class JSONNode:
    def __init__(self, node_id: int, nucleobase: str, column_id: int, block_id: int, aligned_to: int):
        self.id = node_id
        self.nucleobase = nucleobase
        self.column_id = column_id
        self.block_id = block_id
        self.aligned_to = aligned_to


class JSONSequence:
    def __init__(self,
                 node_id: int,
                 genbankID: str,
                 assemblyID: str,
                 mafname: str,
                 name: str,
                 group: str,
                 nodes_ids: List[int]):
        self.id = node_id
        self.genbankID = genbankID
        self.assemblyID = assemblyID
        self.mafname = mafname
        self.name = name
        self.group = group
        self.nodes_ids = nodes_ids


class JSONConsensus:
    def __init__(self,
                 node_id: int,
                 name: str,
                 parent: int,
                 children: List[int],
                 comp_to_all_sequences: Dict[str, float],
                 sequences_ids: List[int],
                 nodes_ids: List[int],
                 mincomp: float):
        self.id = node_id
        self.name = name
        self.parent = parent
        self.children = children
        self.comp_to_all_sequences = comp_to_all_sequences
        self.sequences_ids = sequences_ids
        self.nodes_ids = nodes_ids
        self.mincomp = mincomp


class JSONMAFNode:
    def __init__(self,
                 node_id: int,
                 orient: int,
                 out_edges
                 ):
        self.id = node_id
        self.orient = orient
        self.out_edges = out_edges


class JSONPangenome:
    program_parameters: JSONProgramParameters
    consensuses_tree: List[JSONConsensus]
    sequences: List[JSONSequence]
    nodes: List[JSONNode]

    def __init__(self, pangenome: Pangenome = None, program_parameters: PangenomeParameters = None):
        if program_parameters:
            self.program_parameters = JSONProgramParameters(program_parameters)
        else:
            self.program_parameters = None

        if pangenome.dagmaf:
            self.dagmaf = [JSONMAFNode(node_id=n.id,
                                       orient=n.orient,
                                       out_edges=n.out_edges)
                           for n in pangenome.dagmaf.dagmafnodes]
        else:
            self.dagmaf = []

        if pangenome.pangraph.nodes:
            self.nodes = [JSONNode(node_id=node.id,
                                   nucleobase=node.base.decode("ASCII"),
                                   column_id=node.column_id,
                                   block_id=node.block_id,
                                   aligned_to=node.aligned_to)
                          for node in pangenome.pangraph.nodes]
        else:
            self.nodes = None

        paths_str_id_to_int_id = {seq_id: i for i, seq_id in enumerate(sorted(pangenome.pangraph.paths.keys()))}
        if pangenome.pangraph.paths:
            seqeuences_metadata = [pangenome.genomes_info.genomes_metadata[seqID]
                                   for seqID in pangenome.pangraph.paths.keys()]
            sorted_seqeuences_metadata = sorted(seqeuences_metadata, key=lambda m: paths_str_id_to_int_id[m.mafname])

            self.sequences = [JSONSequence(node_id=paths_str_id_to_int_id[seq_metadata.mafname],  # todo główne ID!!!
                                           genbankID=seq_metadata.genbankID,
                                           assemblyID=seq_metadata.assemblyID,
                                           mafname=seq_metadata.mafname,
                                           name=seq_metadata.name,
                                           group=seq_metadata.group,
                                           nodes_ids=list(itertools.chain.from_iterable(
                                               pangenome.pangraph.paths[seq_metadata.mafname])))
                              for i, seq_metadata in enumerate(sorted_seqeuences_metadata)]
        else:
            self.sequences = None

        if pangenome.consensuses_tree:
            self.consensuses = [JSONConsensus(node_id=consensus_node.consensus_id,
                                              name=f"CONSENSUS{consensus_node.consensus_id}",
                                              parent=consensus_node.parent_node_id,
                                              children=consensus_node.children_nodes_ids,
                                              comp_to_all_sequences=consensus_node.compatibilities_to_all,
                                              sequences_ids=[paths_str_id_to_int_id[seq_id]
                                                             for seq_id in consensus_node.sequences_ids],
                                              nodes_ids=consensus_node.consensus_path,
                                              mincomp=consensus_node.mincomp)
                                for consensus_node in pangenome.consensuses_tree.nodes]
        else:
            self.consensuses = []


