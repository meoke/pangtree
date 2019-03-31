import itertools
from typing import List, Dict
from Pangenome import Pangenome
from PangenomeParameters import PangenomeParameters, FastaComplementationOption
from data.custom_types import NodeID


class JSONProgramParameters:
    def __init__(self, params: PangenomeParameters):
        # input data
        self.multialignment_file_path: str = str(params.multialignment_file_path)
        self.multialignment_format: str = str(params.multialignment_format)

        self.datatype: str = str(params.datatype)
        self.metadata_file_path: str = str(params.metadata_file_path)
        self.blosum_file_path: str = str(params.blosum_file_path)

        # output spec
        self.output_path: str = str(params.output_path)
        self.output_po: bool = params.output_po
        self.generate_fasta: bool = params.generate_fasta
        self.output_with_nodes = params.output_with_nodes
        self.verbose = params.verbose
        self.quiet = params.quiet

        # build spec
        self.raw_maf: bool = params.raw_maf
        self.fasta_complementation_option: FastaComplementationOption = str(params.fasta_complementation_option)
        self.email: str = params.email_address
        self.cache: bool = params.cache
        self.missing_base_symbol: str = params.missing_base_symbol
        self.fasta_source_file: str = params.fasta_source_file

        # consensus spec
        self.consensus_type: str = str(params.consensus_type)
        self.hbmin: float = params.hbmin

        self.max_cutoff_strategy = params.max_cutoff_option
        self.search_range = params.search_range

        self.node_cutoff_strategy = params.node_cutoff_option
        self.multiplier: float = params.multiplier

        self.stop: float = params.stop
        self.re_consensus: bool = params.re_consensus
        self.p: float = params.p


class JSONNode:
    def __init__(self, node_id: int, base: str, column_id: int, block_id: int, aligned_to: int):
        self.id = node_id
        self.base = base
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

        if pangenome.pangraph.nodes and program_parameters.output_with_nodes:
            self.nodes = [JSONNode(node_id=node.id,
                                   base=node.base.decode("ASCII"),
                                   column_id=node.column_id,
                                   block_id=node.block_id,
                                   aligned_to=node.aligned_to)
                          for node in pangenome.pangraph.nodes]
        else:
            self.nodes = None

        paths_str_id_to_int_id = {seq_id: i for i, seq_id in enumerate(sorted(pangenome.pangraph.paths.keys()))}
        if pangenome.pangraph.paths:
            self.sequences = [JSONSequence(sequence_int_id=i,
                                           sequence_str_id=seq_id,
                                           metadata=pangenome.genomes_info.get_seq_metadata_as_dict(seq_id),
                                           nodes_ids=list(itertools.chain.from_iterable(
                                               pangenome.pangraph.paths[seq_id])) if program_parameters.output_with_nodes else []
                                           )
                              for i, seq_id in enumerate(sorted(pangenome.genomes_info.get_all_sequences_ids()))]
        else:
            self.sequences = None

        if pangenome.consensuses_tree:
            self.consensuses = [JSONConsensus(node_id=consensus_node.consensus_id,
                                              name=f"CONSENSUS{consensus_node.consensus_id}",
                                              parent=consensus_node.parent_node_id,
                                              children=consensus_node.children_nodes_ids,
                                              comp_to_all_sequences={seq_id: comp.value
                                                                     for seq_id, comp in
                                                                     consensus_node.compatibilities_to_all.items()},
                                              sequences_ids=[paths_str_id_to_int_id[seq_id]
                                                             for seq_id in
                                                             consensus_node.sequences_ids],
                                              nodes_ids=consensus_node.consensus_path if program_parameters.output_with_nodes else [],
                                              mincomp=consensus_node.mincomp.value)
                                for consensus_node in pangenome.consensuses_tree.nodes]
        else:
            self.consensuses = []


