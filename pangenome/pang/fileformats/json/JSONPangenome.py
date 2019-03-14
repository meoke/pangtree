import itertools
from typing import List, Dict
from Pangenome import Pangenome
from arguments.PangenomeParameters import PangenomeParameters, FastaComplementationOption
from pangraph.custom_types import NodeID


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
                 sequence_int_id: int,
                 sequence_str_id: str,
                 metadata: Dict,
                 nodes_ids: List[NodeID]):
        self.sequence_int_id: int = sequence_int_id
        self.sequence_str_id: str = sequence_str_id
        self.metadata: Dict = metadata
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
        self.node_id = node_id
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
            self.sequences = [JSONSequence(sequence_int_id = i,
                                           sequence_str_id = seq_id,
                                           metadata = pangenome.genomes_info.get_seq_metadata_as_dict(seq_id),
                                           nodes_ids=list(itertools.chain.from_iterable(
                                               pangenome.pangraph.paths[seq_id]))
                                           )
                              for i, seq_id in enumerate(sorted(pangenome.genomes_info.get_all_sequences_ids()))]
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


