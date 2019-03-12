from io import StringIO
from collections import namedtuple
from typing import Dict, Tuple, List, Optional

from mafgraph.graph import Block
from mafgraph.graph.Arc import Arc

from fileformats.maf.DAGMaf import DAGMaf
from fileformats.maf.reader import maf_to_dagmaf
from pangraph.FastaSource import FastaSource
from pangraph.Node import Node
from pangraph.PangraphBuilder.PangraphBuilderBase import PangraphBuilderBase
from pangraph.exceptions import NoSequenceInfo, SequenceBuildingException
from pangraph.custom_types import ColumnID, SequenceID, NodeID, BlockID, Sequence, Nucleobase, make_nucleobase
from metadata.MultialignmentMetadata import MultialignmentMetadata
from mafgraph.mafreader import start_position

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


class PangraphBuilderFromDAG(PangraphBuilderBase):
    def __init__(self,
                 genomes_info: MultialignmentMetadata,
                 missing_nucleotide_symbol: str,
                 fasta_source: Optional[FastaSource] = None):
        super().__init__(genomes_info)
        self.pangraph = None
        self.dagmaf: DAGMaf = None
        self.free_edges: Dict[SequenceID, List[Edge]] = None
        self.seqs_info: Dict[SequenceID, SequenceInfo] = None
        self.column_id: ColumnID = None
        self.complement_sequences: bool = True if fasta_source else False  # todo use this info while building
        self.missing_nucleotide_symbol: Nucleobase = make_nucleobase(missing_nucleotide_symbol)

        if self.complement_sequences:
            self.full_sequences: Dict[SequenceID, Sequence] = self.get_sequences(fasta_source)

    def get_sequences(self, fasta_source: FastaSource) -> Dict[SequenceID, Sequence]:
        return {
            SequenceID(seq_id): Sequence(fasta_source.get_source(sequenceID=seq_id))
            for seq_id in self.sequences_names
        }

    def build(self, maf: StringIO, pangraph) -> None:
        self.dagmaf = maf_to_dagmaf(maf)
        self.pangraph = pangraph
        self.init_pangraph()
        self.free_edges = {SequenceID(seq_id): [] for seq_id in self.sequences_names}
        self.set_seqs_info()
        self.column_id = ColumnID(-1)

        self.complement_starting_nodes()
        for mafnode in self.dagmaf.dagmafnodes:
            self.process_block(mafnode)

    def init_pangraph(self):
        self.pangraph.paths = {seq_id: [] for seq_id in self.sequences_names}

    def set_seqs_info(self) -> None:
        self.seqs_info = {SequenceID(seq_id): [] for seq_id in self.sequences_names}

        for n in self.dagmaf.dagmafnodes:
            for seq in n.alignment:
                self.seqs_info[seq.id].append(SequenceInfo(block_id=BlockID(n.id),
                                                           start=start_position(seq),
                                                           strand=seq.annotations["strand"],
                                                           size=seq.annotations["size"],
                                                           srcSize=seq.annotations["srcSize"],
                                                           orient=n.orient))
        absents_sequences: List[SequenceID] = []
        for seq_id, seq_info_list in self.seqs_info.items():
            if seq_info_list:
                self.seqs_info[seq_id] = sorted(seq_info_list, key=lambda si: si.start)
            else:
                absents_sequences.append(seq_id)
        for seq_id in absents_sequences:
            del self.seqs_info[seq_id]

    def complement_starting_nodes(self) -> None:
        for seq_id, seq_info_list in self.seqs_info.items():
            first_block_sinfo = seq_info_list[0]
            if first_block_sinfo.start != 0:
                self.complement_sequence_starting_nodes(seq_id, first_block_sinfo)

    def get_missing_nucleotide(self, seq_id, i) -> Nucleobase:
        if self.complement_sequences:
            return make_nucleobase(self.full_sequences[seq_id][i])
        return self.missing_nucleotide_symbol

    def complement_sequence_starting_nodes(self, seq_id: SequenceID, first_block_sinfo: SequenceInfo) -> None:
        current_node_id: NodeID = self.get_max_node_id()
        column_id = -first_block_sinfo.start
        join_with = None
        for i in range(first_block_sinfo.start):
            current_node_id += 1
            missing_nucleotide = self.get_missing_nucleotide(seq_id, i)
            self.add_node(node_id=current_node_id,
                          nucleobase=missing_nucleotide,
                          aligned_to=None,
                          column_id=column_id,
                          block_id=None)
            self.add_node_to_sequence(seq_id=seq_id, join_with=join_with, node_id=current_node_id)
            join_with = current_node_id
            column_id += 1
        self.free_edges[seq_id].append(Edge(seq_id=seq_id,
                                            from_block_id=None,
                                            to_block_id=first_block_sinfo.block_id,
                                            last_node_id=current_node_id,))

    def get_max_node_id(self) -> NodeID:
        return NodeID(len(self.pangraph.nodes) - 1)

    def add_node(self,
                 node_id: NodeID,
                 nucleobase: Nucleobase,
                 aligned_to: NodeID,
                 column_id: ColumnID,
                 block_id: BlockID) -> None:
        self.pangraph.nodes.append(Node(node_id=node_id,
                                        base=nucleobase,
                                        aligned_to=aligned_to,
                                        column_id=column_id,
                                        block_id=block_id))

    def add_node_to_sequence(self, seq_id: SequenceID, join_with: NodeID, node_id: NodeID) -> None:
        if not self.pangraph.paths[seq_id] or join_with is None:
            self.pangraph.paths[seq_id].append([node_id])
        else:
            for path in self.pangraph.paths[seq_id]:
                if path[-1] == join_with:
                    path.append(node_id)
                    return
            raise SequenceBuildingException("Cannot find path with specified last node id.")

    def process_block(self, block: Block) -> None:
        print(f"Process block {block.id}")
        current_node_id = self.get_max_node_id()
        block_width = len(block.alignment[0].seq)
        paths_join_info = self.get_paths_join_info(block)

        self.column_id = self.get_max_column_id()
        for col in range(block_width):
            self.column_id += 1
            sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block.alignment}
            nodes_codes = self.get_column_nucleotides_sorted_codes(sequence_name_to_nucleotide)
            column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]
            for i, nucl in enumerate(nodes_codes):
                current_node_id += 1
                seqs_id = [seq_id for seq_id, n in sequence_name_to_nucleotide.items() if n == nucl]
                self.add_node(node_id=current_node_id,
                              nucleobase=make_nucleobase(nucl),
                              aligned_to=self.get_next_aligned_node_id(i, column_nodes_ids),
                              column_id=self.column_id,
                              block_id=block.id)
                for seq_id in seqs_id:
                    self.add_node_to_sequence(seq_id, paths_join_info[seq_id], current_node_id)
                    paths_join_info[seq_id] = current_node_id

        self.add_block_out_edges_to_free_edges(block, paths_join_info)
        self.manage_endings(block, paths_join_info)

    def get_paths_join_info(self, block: Block) -> Dict[SequenceID, NodeID]:
        paths_join_info: Dict[SequenceID, NodeID] = dict()
        for seq in block.alignment:
            paths_join_info[seq.id] = None
            for i, edge in enumerate(self.free_edges[seq.id]):
                if edge.to_block_id == block.id:
                    paths_join_info[seq.id] = edge.last_node_id
        return paths_join_info

    def get_max_column_id(self) -> ColumnID:
        current_columns_ids = [node.column_id for node in self.pangraph.nodes]
        return max(current_columns_ids) if current_columns_ids else ColumnID(-1)

    def get_correct_edge_type(self, block, edge: Arc) -> Tuple[int, int]:
        return edge.edge_type


    def should_join_with_last_node(self, edge_type: Tuple[int, int]) -> bool:
        if edge_type == (1, 1) or edge_type == (1, -1):
            return True
        elif edge_type == (-1, 1) or edge_type == (-1, -1):
            return False
        else:
            raise SequenceBuildingException("Incorrect edge type."
                                            "Cannot decide if sequence should be joined with complemented nucleotides.")

    def complement_sequence_middles_if_needed(self, block: Block, edge: Arc, seq, last_node_id: NodeID):
        seq_id = seq[0].seq_id
        left_block_sinfo, right_block_sinfo = self.get_edge_sinfos(from_block_id=block.id,
                                                                   edge=edge,
                                                                   seq_id=seq_id)
        if self.complementation_not_needed(left_block_sinfo, right_block_sinfo):
            if edge.edge_type == (1, -1):
                return last_node_id
            else:
                return None
        else:
            current_node_id = self.get_max_node_id()
            column_id = self.column_id
            if left_block_sinfo.start < right_block_sinfo.start:
                last_pos = left_block_sinfo.start + left_block_sinfo.size - 1
                next_pos = right_block_sinfo.start
            else:
                last_pos = right_block_sinfo.start + right_block_sinfo.size - 1
                next_pos = left_block_sinfo.start

            join_with = last_node_id if self.should_join_with_last_node(edge.edge_type) else None
            for i in range(last_pos + 1, next_pos):
                column_id += 1
                current_node_id += 1
                missing_nucleotide = self.get_missing_nucleotide(seq_id, i)
                self.add_node(node_id=current_node_id,
                              nucleobase=missing_nucleotide,
                              aligned_to=None,
                              column_id=column_id,
                              block_id=None)
                self.add_node_to_sequence(seq_id=seq_id, join_with=join_with, node_id=current_node_id)
                join_with = current_node_id

            if self.should_join_with_next_node(edge.edge_type):
                return current_node_id
            else:
                return None

    def should_join_with_next_node(self, edge_type):
        if edge_type == (-1, 1) or edge_type == (1, -1) or edge_type == (-1, -1):
            return True
        elif edge_type == (1, 1):
            return False
        else:
            raise SequenceBuildingException("Incorrect edge type."
                                            "Cannot decide if complemented nucleotides must be joined with next block.")

    def add_block_out_edges_to_free_edges(self, block: Block, join_info: Dict[SequenceID, NodeID]):
        for edge in block.out_edges:
            edge_type = self.get_correct_edge_type(block, edge)
            for seq in edge.sequences:
                seq_id = seq[0].seq_id
                last_node_id = self.complement_sequence_middles_if_needed(block=block,
                                                                          edge=edge,
                                                                          seq=seq,
                                                                          last_node_id=join_info[seq_id])
                if last_node_id is not None:
                    self.free_edges[seq_id].append(Edge(seq_id=seq_id,
                                                        from_block_id=block.id,
                                                        to_block_id=edge.to,
                                                        last_node_id=last_node_id))

    def manage_endings(self, block: Block, join_info: Dict[SequenceID, NodeID]):
        sequences_ending_in_this_block = self.get_ending_sequences(block)
        for seq_id in sequences_ending_in_this_block:
            block_sinfo: SequenceInfo = self.get_sinfo(seq_id, block.id)
            if self.sequence_not_complete(block_sinfo):
                last_node_id = self.complement_sequence_middle_nodes(seq_id=seq_id,
                                                                     last_pos=block_sinfo.start + block_sinfo.size-1,
                                                                     next_pos=block_sinfo.srcSize,
                                                                     last_node_id=join_info[seq_id])
            else:
                last_node_id = join_info[seq_id]
            self.free_edges[seq_id].append(Edge(seq_id=seq_id,
                                                from_block_id=block.id,
                                                to_block_id=None,
                                                last_node_id=last_node_id))

    def get_edge_sinfos(self, from_block_id: BlockID, edge: Arc, seq_id: SequenceID) -> \
            Tuple[SequenceInfo, SequenceInfo]:
        left_seq_info, right_seq_info = None, None
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == from_block_id:
                left_seq_info = sinfo
            if sinfo.block_id == edge.to:
                right_seq_info = sinfo
        if left_seq_info is None or right_seq_info is None:
            raise NoSequenceInfo(f"SequenceInfos for edge cannot be None. "
                                 f"Left block is {left_seq_info}, right block is {right_seq_info}.")
        return left_seq_info, right_seq_info

    def get_column_nucleotides_sorted_codes(self, seq_to_nucl: Dict[SequenceID, str]) -> List[str]:
        return sorted(
                set(
                    [nucleotide
                     for nucleotide
                     in seq_to_nucl.values()
                     if nucleotide is not '-']))

    def get_sinfo(self, seq_id, block_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == block_id:
                return sinfo

    def complementation_not_needed(self, left, right):
        return left.start + left.size == right.start or right.start + right.size == left.start

    def complement_sequence_middle_nodes(self, seq_id: SequenceID, last_pos, next_pos, last_node_id: NodeID) -> NodeID:
        current_node_id = self.get_max_node_id()
        column_id = self.column_id
        join_with = last_node_id
        for i in range(last_pos+1, next_pos):
            column_id += 1
            current_node_id += 1
            missing_nucleotide = self.get_missing_nucleotide(seq_id, i)
            self.add_node(node_id=current_node_id,
                          nucleobase=missing_nucleotide,
                          aligned_to=None,
                          column_id=column_id,
                          block_id=None)
            self.add_node_to_sequence(seq_id=seq_id, join_with=join_with, node_id=current_node_id)
            join_with = current_node_id
        return current_node_id

    def get_ending_sequences(self, block):
        sequences_ending_in_this_block = []
        for seq_id, sinfo_list in self.seqs_info.items():
            last_block_sinfo = sinfo_list[-1]
            if last_block_sinfo.block_id == block.id:
                sequences_ending_in_this_block.append(seq_id)
        return sequences_ending_in_this_block

    def sequence_not_complete(self, last_block_sinfo):
        if last_block_sinfo.strand == 1:
            return last_block_sinfo.start + last_block_sinfo.size != last_block_sinfo.srcSize
        elif last_block_sinfo.strand == -1:
            return last_block_sinfo.start != 0
        else:
            raise Exception("Unecpected strand value")

    def get_next_aligned_node_id(self, current_column_i, column_nodes_ids):
        if len(column_nodes_ids) > 1:
            return column_nodes_ids[(current_column_i + 1) % len(column_nodes_ids)]
        return None

