from collections import namedtuple
from typing import Tuple, List, NewType, Optional, Dict
from datamodel.DAGMaf import DAGMaf
from datamodel.Sequence import Sequences, SequenceID, Sequence, SequencePath
from datamodel.Node import Node, ColumnID, BlockID, NodeID, Base
from datamodel.builders.PoagraphBuildException import PoagraphBuildException
from datamodel.fasta_providers import FastaProvider
from datamodel.input_types import MetadataCSV

from mafgraph.mafreader import start_position

MafSequenceID = NewType('MafSequenceID', str)

SequenceInfo = namedtuple('SequenceInfo', ['block_id',
                                           'start',
                                           'strand',
                                           'size',
                                           'srcSize',
                                           'orient'])

Edge = namedtuple('Edge', ['seq_id',
                           'from_block_id',
                           'to_block_id',
                           'last_node_id'])

class _BuildState:
    def __init__(self,
                 initial_nodes: List[Node],
                 initial_sequences: Sequences,
                 initial_edges: Dict[SequenceID, List[Edge]],
                 seqs_info: Dict[SequenceID, SequenceInfo],
                 initial_column_id: ColumnID,
                 fasta_provider: FastaProvider):
        self.nodes: List[Node] = initial_nodes
        self.sequences: Sequences = initial_sequences
        self.free_edges: Dict[SequenceID, List[Edge]] = initial_edges
        self.seqs_info: Dict[SequenceID, SequenceInfo] = seqs_info
        self.column_id: ColumnID = initial_column_id
        self.fasta_provider: FastaProvider = fasta_provider


def get_poagraph(dagmaf: DAGMaf,
                 fasta_provider: FastaProvider,
                 metadata: Optional[MetadataCSV]) -> Tuple[List[Node], Sequences]:
    sequences_in_dagmaf = _get_sequences_ids(dagmaf)
    build_state = _BuildState(initial_nodes=[],
                              initial_sequences=_init_sequences(sequences_in_dagmaf, metadata),
                              initial_edges=_init_free_edges(sequences_in_dagmaf),
                              seqs_info=_get_seqs_info(dagmaf, sequences_in_dagmaf),
                              initial_column_id=ColumnID(-1),
                              fasta_provider=fasta_provider)

    _complement_starting_nodes(build_state)


def _get_sequences_ids(dagmaf: DAGMaf) -> List[SequenceID]:
    return list({SequenceID(seq.id) for block in dagmaf.dagmafnodes for seq in block.alignment})


def _init_sequences(sequences_in_dagmaf: List[SequenceID],
                    metadata: Optional[MetadataCSV]) -> Sequences:
    metadata_sequences_ids = metadata.get_all_sequences_ids() if metadata else []
    initial_sequences: Sequences = Sequences({seq_id: Sequence(seqid=seq_id,
                                                               paths=[],
                                                               seqmetadata=metadata.get_sequence_metadata(seq_id)
                                                               if metadata else {})
                                              for seq_id in set(sequences_in_dagmaf + metadata_sequences_ids)})

    return initial_sequences


def _init_free_edges(maf_sequences_ids: List[SequenceID]) -> Dict[SequenceID, List[Edge]]:
    return {seq_id: [] for seq_id in maf_sequences_ids}


def _get_seqs_info(dagmaf: DAGMaf, sequences_in_dagmaf: List[SequenceID]) -> Dict[SequenceID, SequenceInfo]:
    seqs_info = {seq_id: [] for seq_id in sequences_in_dagmaf}

    for n in dagmaf.dagmafnodes:
        for seq in n.alignment:
            seqs_info[SequenceID(seq.id)].append(SequenceInfo(block_id=BlockID(n.id),
                                                              start=start_position(seq),
                                                              strand=seq.annotations["strand"],
                                                              size=seq.annotations["size"],
                                                              srcSize=seq.annotations["srcSize"],
                                                              orient=n.orient))
    absents_sequences: List[SequenceID] = []
    for seq_id, seq_info_list in seqs_info.items():
        if seq_info_list:
            seqs_info[seq_id] = sorted(seq_info_list, key=lambda si: si.start)
        else:
            absents_sequences.append(seq_id)
    for seq_id in absents_sequences:
        del seqs_info[seq_id]
    return seqs_info


def _complement_starting_nodes(build_state: _BuildState) -> None:
    for seq_id, seq_info_list in build_state.seqs_info.items():
        first_block_sinfo = seq_info_list[0]
        if first_block_sinfo.start != 0:
            _complement_sequence_starting_nodes(build_state, seq_id, first_block_sinfo)


def _complement_sequence_starting_nodes(build_state: _BuildState,
                                        seq_id: SequenceID,
                                        first_block_sinfo: SequenceInfo) -> None:
    current_node_id: NodeID = _get_max_node_id(build_state.nodes)
    column_id = -first_block_sinfo.start
    join_with = None
    for i in range(first_block_sinfo.start):
        current_node_id += 1
        missing_nucleotide = _get_missing_nucleotide(seq_id, i, build_state.fasta_provider)
        build_state.nodes += Node(node_id=current_node_id,
                                  base=missing_nucleotide,
                                  aligned_to=None,
                                  column_id=column_id)
        _add_node_to_sequence(build_state, seq_id=seq_id, join_with=join_with, node_id=current_node_id)
        join_with = current_node_id
        column_id += 1
    build_state.free_edges + Edge(seq_id=seq_id,
                                  from_block_id=None,
                                  to_block_id=first_block_sinfo.block_id,
                                  last_node_id=current_node_id)


def _get_max_node_id(nodes: List[Node]) -> NodeID:
    return NodeID(len(nodes) - 1)


def _get_missing_nucleotide(fasta_provider, seq_id: SequenceID, i: int) -> Base:
    return fasta_provider.get_base(seq_id, i)


def _add_node_to_sequence(build_state: _BuildState,
                          seq_id: SequenceID,
                          join_with: NodeID,
                          node_id: NodeID) -> None:
    if not build_state.sequences[seq_id].paths or join_with is None:
        build_state.sequences[seq_id].paths.append(SequencePath([node_id]))
    else:
        for path in build_state.sequences[seq_id].paths:
            if path[-1] == join_with:
                path.append(node_id)
                return
    raise PoagraphBuildException("Cannot find path with specified last node id.")