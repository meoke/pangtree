from collections import namedtuple
from typing import Tuple, List, NewType, Optional, Dict

from pangtreebuild.mafgraph.graph import Block
from pangtreebuild.mafgraph.graph.Arc import Arc
from pangtreebuild.mafgraph.mafreader import start_position

from pangtreebuild.poagraph.DAGMaf import DAGMaf, DAGMafNode
from pangtreebuild.poagraph.Sequence import SequenceID, Sequence, SeqPath
from pangtreebuild.poagraph.Node import Node, ColumnID, BlockID, NodeID, Base
from pangtreebuild.poagraph.fasta_providers import FastaProvider
from pangtreebuild.poagraph.input_types import MetadataCSV

from pangtreebuild.tools import logprocess



global_logger = logprocess.get_global_logger()
detailed_logger = logprocess.get_logger("details")


class PoagraphBuildException(Exception):
    pass


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
                 initial_sequences: Dict[SequenceID, Sequence],
                 initial_edges: Dict[SequenceID, List[Edge]],
                 seqs_info: Dict[SequenceID, List[SequenceInfo]],
                 initial_column_id: ColumnID,
                 fasta_provider: FastaProvider):
        self.nodes: List[Node] = initial_nodes
        self.sequences: Dict[SequenceID, Sequence] = initial_sequences
        self.free_edges: Dict[SequenceID, List[Edge]] = initial_edges
        self.seqs_info: Dict[SequenceID, List[SequenceInfo]] = seqs_info
        self.column_id: ColumnID = initial_column_id
        self.fasta_provider: FastaProvider = fasta_provider


def get_poagraph(dagmaf: DAGMaf,
                 fasta_provider: FastaProvider,
                 metadata: Optional[MetadataCSV]) -> Tuple[List[Node], Dict[SequenceID, Sequence]]:
    sequences_in_dagmaf = _get_sequences_ids(dagmaf)
    build_state = _BuildState(initial_nodes=[],
                              initial_sequences=_init_sequences(sequences_in_dagmaf, metadata),
                              initial_edges=_init_free_edges(sequences_in_dagmaf),
                              seqs_info=_get_seqs_info(dagmaf, sequences_in_dagmaf),
                              initial_column_id=ColumnID(-1),
                              fasta_provider=fasta_provider)

    _complement_starting_nodes(build_state)

    for i, mafnode in enumerate(dagmaf.dagmaf_nodes):
        # p = i/len(dagmaf.dagmaf_nodes) * 100
        # global_logger.info(f"{p}%")
        _process_block(build_state, mafnode)

    return build_state.nodes, build_state.sequences


def _get_sequences_ids(dagmaf: DAGMaf) -> List[SequenceID]:
    return list({SequenceID(seq.id) for block in dagmaf.dagmaf_nodes for seq in block.alignment})


def _init_sequences(sequences_in_dagmaf: List[SequenceID],
                    metadata: Optional[MetadataCSV]) -> Dict[SequenceID, Sequence]:
    metadata_sequences_ids = metadata.get_all_sequences_ids() if metadata else []
    initial_sequences = {seq_id: Sequence(seqid=seq_id,
                                          paths=[],
                                          seqmetadata=metadata.get_sequence_metadata(seq_id)
                                          if metadata else {})
                         for seq_id in set(sequences_in_dagmaf + metadata_sequences_ids)}

    return initial_sequences


def _init_free_edges(maf_sequences_ids: List[SequenceID]) -> Dict[SequenceID, List[Edge]]:
    return {seq_id: [] for seq_id in maf_sequences_ids}


def _get_seqs_info(dagmaf: DAGMaf, sequences_in_dagmaf: List[SequenceID]) -> Dict[SequenceID, List[SequenceInfo]]:
    seqs_info = {seq_id: [] for seq_id in sequences_in_dagmaf}

    for n in dagmaf.dagmaf_nodes:
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
        missing_nucleotide = _get_missing_nucleotide(build_state.fasta_provider, seq_id, i)
        build_state.nodes += [Node(node_id=current_node_id,
                                   base=missing_nucleotide,
                                   column_id=column_id)]
        _add_node_to_sequence(build_state, seq_id=seq_id, join_with=join_with, node_id=current_node_id)
        join_with = current_node_id
        column_id += 1
    build_state.free_edges[seq_id] += [Edge(seq_id=seq_id,
                                            from_block_id=None,
                                            to_block_id=first_block_sinfo.block_id,
                                            last_node_id=current_node_id)]


def _get_max_node_id(nodes: List[Node]) -> NodeID:
    return NodeID(len(nodes) - 1)


def _get_missing_nucleotide(fasta_provider, seq_id: SequenceID, i: int) -> Base:
    return fasta_provider.get_base(seq_id, i)


def _add_node_to_sequence(build_state: _BuildState,
                          seq_id: SequenceID,
                          join_with: NodeID,
                          node_id: NodeID) -> None:
    q = build_state.sequences[seq_id].paths
    if len(build_state.sequences[seq_id].paths) == 0 or join_with is None:
        build_state.sequences[seq_id].paths.append(SeqPath([node_id]))
    else:
        for path in build_state.sequences[seq_id].paths:
            if path[-1] == join_with:
                path.append(node_id)
                return

        raise PoagraphBuildException("Cannot find path with specified last node id_.")


def _process_block(build_state: _BuildState, block: DAGMafNode):
    # global_logger.info(f"Processing block {block.id_}...")
    current_node_id = _get_max_node_id(build_state.nodes)
    block_width = len(block.alignment[0].seq)
    paths_join_info = _get_paths_join_info(block, build_state.free_edges)

    build_state.column_id = _get_max_column_id(build_state.nodes)
    for col in range(block_width):
        build_state.column_id += 1
        sequence_name_to_nucleotide = {MafSequenceID(seq.id): seq[col] for seq in block.alignment}
        nodes_codes = _get_column_nucleotides_sorted_codes(sequence_name_to_nucleotide)
        column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]
        for i, nucl in enumerate(nodes_codes):
            current_node_id += 1
            maf_seqs_id = [seq_id for seq_id, n in sequence_name_to_nucleotide.items() if n == nucl]
            build_state.nodes += [Node(node_id=current_node_id,
                                       base=Base(nucl),
                                       aligned_to=_get_next_aligned_node_id(i, column_nodes_ids),
                                       column_id=build_state.column_id,
                                       block_id=block.id)]

            for maf_seq_id in maf_seqs_id:
                seq_id = SequenceID(maf_seq_id)
                _add_node_to_sequence(build_state, seq_id, paths_join_info[seq_id], current_node_id)
                paths_join_info[seq_id] = current_node_id

    _add_block_out_edges_to_free_edges(build_state, block, paths_join_info)
    _manage_endings(build_state, block, paths_join_info)


def _get_paths_join_info(block: Block,
                         free_edges: Dict[SequenceID, List[Edge]]) -> Dict[SequenceID, Optional[NodeID]]:
    paths_join_info: Dict[SequenceID, Optional[NodeID]] = dict()
    for seq in block.alignment:
        seq_id = SequenceID(seq.id)
        paths_join_info[seq_id] = None
        for i, edge in enumerate(free_edges[seq_id]):
            if edge.to_block_id == block.id:
                paths_join_info[seq_id] = edge.last_node_id
    return paths_join_info


def _get_max_column_id(nodes: List[Node]) -> ColumnID:
    current_columns_ids = [node.column_id for node in nodes]
    return max(current_columns_ids) if current_columns_ids else ColumnID(-1)


def _get_column_nucleotides_sorted_codes(seq_to_nucl: Dict[SequenceID, str]) -> List[str]:
    return sorted(
        set(
            [nucleotide
             for nucleotide
             in seq_to_nucl.values()
             if nucleotide is not '-']))


def _get_next_aligned_node_id(current_column_i, column_nodes_ids) -> Optional[NodeID]:
    if len(column_nodes_ids) > 1:
        return column_nodes_ids[(current_column_i + 1) % len(column_nodes_ids)]
    return None


def _add_block_out_edges_to_free_edges(build_state: _BuildState,
                                       block: Block,
                                       join_info: Dict[SequenceID, NodeID]):
    for edge in block.out_edges:
        _ = _get_correct_edge_type(edge)
        for seq in edge.sequences:
            seq_id = SequenceID(seq[0].seq_id)
            last_node_id = _complement_sequence_middles_if_needed(build_state=build_state,
                                                                  block=block,
                                                                  edge=edge,
                                                                  seq=seq,
                                                                  last_node_id=join_info[seq_id])
            if last_node_id is not None:
                build_state.free_edges[seq_id].append(Edge(seq_id=seq_id,
                                                           from_block_id=block.id,
                                                           to_block_id=edge.to,
                                                           last_node_id=last_node_id))


def _get_correct_edge_type(edge: Arc) -> Tuple[int, int]:
    return edge.edge_type


def _complement_sequence_middles_if_needed(build_state: _BuildState,
                                           block: Block,
                                           edge: Arc,
                                           seq,
                                           last_node_id: NodeID):
    seq_id = SequenceID(seq[0].seq_id)
    left_block_sinfo, right_block_sinfo = _get_edge_sinfos(seqs_info=build_state.seqs_info,
                                                           from_block_id=block.id,
                                                           edge=edge,
                                                           seq_id=seq_id)
    if _complementation_not_needed(left_block_sinfo, right_block_sinfo):
        if edge.edge_type == (1, -1):
            return last_node_id
        else:
            return None
    else:
        current_node_id = _get_max_node_id(build_state.nodes)
        column_id = build_state.column_id
        if left_block_sinfo.start < right_block_sinfo.start:
            last_pos = left_block_sinfo.start + left_block_sinfo.size - 1
            next_pos = right_block_sinfo.start
        else:
            last_pos = right_block_sinfo.start + right_block_sinfo.size - 1
            next_pos = left_block_sinfo.start

        join_with = last_node_id if _should_join_with_last_node(edge.edge_type) else None
        for i in range(last_pos + 1, next_pos):
            column_id += 1
            current_node_id += 1
            missing_nucleotide = _get_missing_nucleotide(build_state.fasta_provider, seq_id, i)
            build_state.nodes += [Node(node_id=current_node_id,
                                       base=missing_nucleotide,
                                       aligned_to=None,
                                       column_id=column_id,
                                       block_id=None)]
            _add_node_to_sequence(build_state, seq_id=seq_id, join_with=join_with, node_id=current_node_id)
            join_with = current_node_id

        if _should_join_with_next_node(edge.edge_type):
            return current_node_id
        else:
            return None


def _get_edge_sinfos(seqs_info: Dict[SequenceID, List[SequenceInfo]],
                     from_block_id: BlockID,
                     edge: Arc,
                     seq_id: SequenceID) -> Tuple[SequenceInfo, SequenceInfo]:
    left_seq_info, right_seq_info = None, None
    for sinfo in seqs_info[seq_id]:
        if sinfo.block_id == from_block_id:
            left_seq_info = sinfo
        if sinfo.block_id == edge.to:
            right_seq_info = sinfo
    if left_seq_info is None or right_seq_info is None:
        raise PoagraphBuildException(f"SequenceInfos for edge cannot be None. "
                                     f"Left block is {left_seq_info}, right block is {right_seq_info}.")
    return left_seq_info, right_seq_info


def _complementation_not_needed(left: SequenceInfo, right: SequenceInfo) -> bool:
    return left.start + left.size == right.start or right.start + right.size == left.start


def _should_join_with_last_node(edge_type: Tuple[int, int]) -> bool:
    if edge_type == (1, 1) or edge_type == (1, -1):
        return True
    elif edge_type == (-1, 1) or edge_type == (-1, -1):
        return False
    else:
        raise PoagraphBuildException("Incorrect edge type."
                                     "Cannot decide if sequence should be joined with complemented nucleotides.")


def _should_join_with_next_node(edge_type: Tuple[int, int]) -> bool:
    if edge_type == (-1, 1) or edge_type == (1, -1) or edge_type == (-1, -1):
        return True
    elif edge_type == (1, 1):
        return False
    else:
        raise PoagraphBuildException("Incorrect edge type."
                                     "Cannot decide if complemented nucleotides must be joined with next block.")


def _manage_endings(build_state: _BuildState, block: Block, join_info: Dict[SequenceID, NodeID]):
    sequences_ending_in_this_block = _get_ending_sequences(build_state.seqs_info, block)
    for seq_id in sequences_ending_in_this_block:
        block_sinfo: SequenceInfo = _get_sinfo(build_state.seqs_info[seq_id], block.id)
        if _sequence_not_complete(block_sinfo):
            last_node_id = _complement_sequence_middle_nodes(build_state,
                                                             seq_id=seq_id,
                                                             last_pos=block_sinfo.start + block_sinfo.size-1,
                                                             next_pos=block_sinfo.srcSize,
                                                             last_node_id=join_info[seq_id])
        else:
            last_node_id = join_info[seq_id]
        build_state.free_edges[seq_id].append(Edge(seq_id=seq_id,
                                                   from_block_id=block.id,
                                                   to_block_id=None,
                                                   last_node_id=last_node_id))


def _get_ending_sequences(seqs_info: Dict[SequenceID, List[SequenceInfo]], block: Block) -> List[SequenceID]:
    sequences_ending_in_this_block = []
    for seq_id, sinfo_list in seqs_info.items():
        last_block_sinfo = sinfo_list[-1]
        if last_block_sinfo.block_id == block.id:
            sequences_ending_in_this_block.append(seq_id)
    return sequences_ending_in_this_block


def _get_sinfo(seq_info: List[SequenceInfo], block_id: int) -> SequenceInfo:
    for sinfo in seq_info:
        if sinfo.block_id == block_id:
            return sinfo
    raise PoagraphBuildException(f"No sequences info for given block")


def _sequence_not_complete(last_block_sinfo: SequenceInfo) -> bool:
    if last_block_sinfo.strand == 1:
        return last_block_sinfo.start + last_block_sinfo.size != last_block_sinfo.srcSize
    elif last_block_sinfo.strand == -1:
        return last_block_sinfo.start != 0
    else:
        raise Exception("Unexpected strand value")


def _complement_sequence_middle_nodes(build_state: _BuildState,
                                      seq_id: SequenceID,
                                      last_pos,
                                      next_pos,
                                      last_node_id: NodeID) -> NodeID:
    current_node_id = _get_max_node_id(build_state.nodes)
    column_id = build_state.column_id
    join_with = last_node_id
    for i in range(last_pos+1, next_pos):
        column_id += 1
        current_node_id += 1
        missing_nucleotide = _get_missing_nucleotide(build_state.fasta_provider, seq_id, i)
        build_state.nodes += [Node(node_id=current_node_id,
                                   base=missing_nucleotide,
                                   aligned_to=None,
                                   column_id=column_id,
                                   block_id=None)
                              ]
        _add_node_to_sequence(build_state, seq_id=seq_id, join_with=join_with, node_id=current_node_id)
        join_with = current_node_id
    return current_node_id
