import sys
from collections import deque, namedtuple

from mafgraph.graph import Block

from fileformats.maf.DAGMaf import DAGMaf
from fileformats.maf.reader import maf_to_dagmaf
from graph import nucleotides
from graph.FastaSource import FastaSource
from graph.Node import Node
from graph.PangraphBuilder.PangraphBuilderBase import PangraphBuilderBase
from metadata.MultialignmentMetadata import MultialignmentMetadata
from mafgraph.mafreader import start_position


class PangraphBuilderFromDAG(PangraphBuilderBase):
    def __init__(self, genomes_info: MultialignmentMetadata, fasta_source: FastaSource):
        super().__init__(genomes_info)
        self.pangraph = None
        self.dagmaf = None
        self.blocks_info = None
        self.free_edges = None
        self.seqs_info = None

        if fasta_source is not None:
            self.full_sequences = self.get_sequences(fasta_source)
            self.complement_sequences = True
        else:
            self.complement_sequences = False

    def get_sequences(self, fasta_source):
        return {
            seq_id: fasta_source.get_source(id=seq_id)
            for seq_id in self.sequences_names
        }

    def get_max_node_id(self):
        current_nodes_ids = [node.id for node in self.pangraph._nodes if node is not None]
        return max(current_nodes_ids) if current_nodes_ids else -1

    def init_pangraph(self):
        nodes_count = PangraphBuilderFromDAG.get_nodes_count(self.dagmaf)
        self.pangraph._nodes = [None] * nodes_count
        self.pangraph._pathmanager.init_paths(self.sequences_names, nodes_count)

    SInfo = namedtuple('SInfo', ['block_id',
                                 'start',
                                 'strand',
                                 'size',
                                 'srcSize'])

    Edge = namedtuple('Edge', ['seq_id', 'from_block_id', 'last_node_id', 'to_block_id'])

    def build(self, input, pangraph):
        self.dagmaf = maf_to_dagmaf(input)
        self.pangraph = pangraph
        self.init_pangraph()

        self.free_edges = {seq_id: [] for seq_id in self.sequences_names}
        self.set_seqs_info()

        self.complement_starting_nodes()
        for mafnode in self.dagmaf.dagmafnodes:
            self.process_block(mafnode)
        # self.complement_ending_nodes()
        self.pangraph.remove_empty_paths()
        self.pangraph.remove_empty_nodes()

    def set_seqs_info(self):
        self.seqs_info = {seq_id: [] for seq_id in self.sequences_names}

        for n in self.dagmaf.dagmafnodes:
            for seq in n.alignment:
                self.seqs_info[seq.id].append(PangraphBuilderFromDAG.SInfo(block_id=n.id,
                                                                           start=start_position(seq),
                                                                           strand=seq.annotations["strand"],
                                                                           size=seq.annotations["size"],
                                                                           srcSize=seq.annotations["srcSize"]))
        to_remove = []
        for seq_id, seq_info_list in self.seqs_info.items():
            if seq_info_list:
                self.seqs_info[seq_id] = sorted(seq_info_list, key=lambda si: si.start)
            else:
                to_remove.append(seq_id)
        for i in to_remove:
            del self.seqs_info[i]

    def get_left_sinfo(self, edge, seq_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == edge.from_block_id:
                return sinfo

    def get_edge_sinfos(self, from_block_id, edge, seq_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == from_block_id:
                left_seq_info = sinfo
            if sinfo.block_id == edge.to:
                right_seq_info = sinfo
        return left_seq_info, right_seq_info

    def complement_starting_nodes(self):
        for seq_id, seq_info_list in self.seqs_info.items():
            first_block_sinfo = seq_info_list[0]
            if first_block_sinfo.start != 0:
                self.complement_sequence_starting_nodes(seq_id, first_block_sinfo)

    def complement_ending_nodes(self):
        for seq_id, edges in self.free_edges.items():
            for edge in edges:
                if edge.to_block_id is None:
                    self.complement_sequence_ending_nodes(seq_id, edge)
                continue

    def complement_sequence_starting_nodes(self, seq_id, first_block_sinfo):
        current_node_id = self.get_max_node_id()
        in_nodes = []
        for i in range(first_block_sinfo.start):
            current_node_id += 1
            missing_nucleotide = self.full_sequences[seq_id][i]
            self.add_node(seqs_id=[seq_id],
                          id=current_node_id,
                          base=missing_nucleotide,
                          in_nodes=in_nodes,
                          aligned_to=None)
            in_nodes = [current_node_id]
        self.free_edges[seq_id].append(PangraphBuilderFromDAG.Edge(seq_id=seq_id,
                                                                   from_block_id=None,
                                                                   last_node_id=current_node_id,
                                                                   to_block_id=first_block_sinfo.block_id))

    def complement_sequence_ending_nodes(self, seq_id, last_edge: Edge):
        current_node_id = self.get_max_node_id()
        in_nodes = [last_edge.last_node_id]
        last_block_sinfo = self.get_sinfo(seq_id, last_edge.from_block_id)
        # if last_block_sinfo.from_block_id == last_edge.from_block_id: # czy może być inaczej?

        for i in range(last_block_sinfo.start + last_block_sinfo.size, last_block_sinfo.srcSize):
            current_node_id += 1
            missing_nucleotide = self.full_sequences[seq_id][i]
            self.add_node(seqs_id=[seq_id],
                          id=current_node_id,
                          base=missing_nucleotide,
                          in_nodes=in_nodes,
                          aligned_to=None)
            in_nodes = [current_node_id]

    def add_node(self, seqs_id, id, base, in_nodes, aligned_to):
        self.pangraph._nodes[id] = Node(id=id,
                                        base=nucleotides.code(base),
                                        in_nodes=in_nodes,
                                        aligned_to=aligned_to)
        for seq_id in seqs_id:
            self.pangraph.add_path_to_node(path_name=seq_id, node_id=id)

    def get_block_in_nodes(self, block):
        block_in_nodes = dict()
        for seq in block.alignment:
            block_in_nodes[seq.id] = None
            for i, edge in enumerate(self.free_edges[seq.id]):
                if edge.to_block_id == block.id:
                    block_in_nodes[seq.id] = edge.last_node_id
                    continue
        return block_in_nodes

    @staticmethod
    def get_column_nucleotides_codes(sequence_name_to_nucleotide):
        return sorted(
                [*(set(
                    [nucleotide
                     for nucleotide
                     in sequence_name_to_nucleotide.values()])).
                    difference(set(['-']))])

    @staticmethod
    def get_in_nodes(block_in_nodes, seqs_id):
        return list(set([node_id for seq_id, node_id in block_in_nodes.items() if seq_id in seqs_id and node_id is not None]))

    @staticmethod
    def update_block_in_nodes(block_in_nodes, seqs_id, current_node_id):
        for seq_id in seqs_id:
            block_in_nodes[seq_id] = current_node_id

    def get_sinfo(self, seq_id, block_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == block_id:
                return sinfo

    @staticmethod
    def continuous_sequence(left_block_sinfo, right_block_sinfo):
        if left_block_sinfo.strand == 1 and right_block_sinfo.strand == 1:
            return left_block_sinfo.start + left_block_sinfo.size == right_block_sinfo.start
        elif left_block_sinfo.strand == 1 and right_block_sinfo.strand == -1:
            raise Exception("+ łączy się z -")
        elif left_block_sinfo.strand == -1 and right_block_sinfo.strand == +1:
            raise Exception("- łączy się z +")
        elif left_block_sinfo.strand == -1 and right_block_sinfo.strand == -1:
            return right_block_sinfo.start + right_block_sinfo.size == left_block_sinfo.start
        else:
            raise Exception("Unexcpected strand values!")

    def complement_sequence_middle_nodes(self, seq_id, last_pos, next_pos, in_node_id) -> int:
        current_node_id = self.get_max_node_id()
        in_nodes = [in_node_id]
        for i in range(last_pos+1, next_pos):
            current_node_id += 1
            missing_nucleotide = self.full_sequences[seq_id][i]
            self.add_node(seqs_id=[seq_id],
                          id=current_node_id,
                          base=missing_nucleotide,
                          in_nodes=in_nodes,
                          aligned_to=None)
            in_nodes = [current_node_id]
        return current_node_id

    def add_block_out_edges_to_free_edges(self, block_in_nodes, block: Block):
        # for edge in block.out_edges:
        #     # if edge.edge_type != (1, -1):
        #     #     continue
        #     for seq in edge.sequences:
        #         seq_id = seq[0].seq_id
        #         left_block_sinfo, right_block_sinfo = self.get_edge_sinfos(from_block_id=block.id, edge=edge, seq_id=seq_id)
        #
        #         if not PangraphBuilderFromDAG.continuous_sequence(left_block_sinfo, right_block_sinfo):
        #             last_node_id = self.complement_sequence_middle_nodes(seq_id=seq_id,
        #                                                                  last_pos=left_block_sinfo.start + left_block_sinfo.size-1,
        #                                                                  next_pos=right_block_sinfo.start,
        #                                                                  in_node_id=block_in_nodes[seq_id])
        #         else:
        #             last_node_id = block_in_nodes[seq_id]
        #             if edge.edge_type != (1, -1):
        #                 continue
        #
        #         if edge.edge_type != (1, -1):
        #             from_block_id = None
        #         else:
        #             from_block_id = block.id
        #
        #         self.free_edges[seq_id].append(
        #             PangraphBuilderFromDAG.Edge(
        #                 seq_id=seq_id,
        #                 from_block_id=from_block_id,
        #                 last_node_id=last_node_id,
        #                 to_block_id=edge.to))

        for edge in block.out_edges:
            if edge.edge_type != (1, -1):
                continue
            for seq in edge.sequences:
                seq_id = seq[0].seq_id
                left_block_sinfo, right_block_sinfo = self.get_edge_sinfos(from_block_id=block.id, edge=edge, seq_id=seq_id)
                if not PangraphBuilderFromDAG.continuous_sequence(left_block_sinfo, right_block_sinfo):
                    last_node_id = self.complement_sequence_middle_nodes(seq_id=seq_id,
                                                                         last_pos=left_block_sinfo.start + left_block_sinfo.size-1,
                                                                         next_pos=right_block_sinfo.start,
                                                                         in_node_id=block_in_nodes[seq_id])
                else:
                    last_node_id = block_in_nodes[seq_id]

                self.free_edges[seq_id].append(
                    PangraphBuilderFromDAG.Edge(
                        seq_id=seq_id,
                        from_block_id=block.id,
                        last_node_id=last_node_id,
                        to_block_id=edge.to))

    def get_ending_sequences(self, block):
        sequences_ending_in_this_block = []
        for seq_id, sinfo_list in self.seqs_info.items():
            last_block_sinfo = sinfo_list[-1]
            if last_block_sinfo.block_id == block.id:
                sequences_ending_in_this_block.append(seq_id)
        return sequences_ending_in_this_block

    def manage_endings(self, block_in_nodes, block):
        sequences_ending_in_this_block = self.get_ending_sequences(block)
        for seq_id in sequences_ending_in_this_block:
            block_sinfo: PangraphBuilderFromDAG.SInfo = self.get_sinfo(seq_id, block.id)
            if self.sequence_not_complete(block_sinfo):
                last_node_id = self.complement_sequence_middle_nodes(seq_id=seq_id,
                                                                     last_pos=block_sinfo.start + block_sinfo.size-1,
                                                                     next_pos=block_sinfo.srcSize,
                                                                     in_node_id=block_in_nodes[seq_id])
            else:
                last_node_id = block_in_nodes[seq_id]
            self.free_edges[seq_id].append(PangraphBuilderFromDAG.Edge(seq_id=seq_id,
                                                                       from_block_id=block.id,
                                                                       last_node_id=last_node_id,
                                                                       to_block_id=None))

    def sequence_not_complete(self, last_block_sinfo):
        if last_block_sinfo.strand == 1:
            return last_block_sinfo.start + last_block_sinfo.size != last_block_sinfo.srcSize
        elif last_block_sinfo.strand == -1:
            return last_block_sinfo.start != 0
        else:
            raise Exception("Unecpected strand value")

    def process_block(self, block):
        current_node_id = self.get_max_node_id()
        block_width = len(block.alignment[0].seq)
        block_in_nodes = self.get_block_in_nodes(block)

        for col in range(block_width):
            sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block.alignment}
            nodes_codes = PangraphBuilderFromDAG.get_column_nucleotides_codes(sequence_name_to_nucleotide)
            column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]
            for i, nucl in enumerate(nodes_codes):
                current_node_id += 1
                seqs_id = [seq_id for seq_id, n in sequence_name_to_nucleotide.items() if n == nucl]
                in_nodes = self.get_in_nodes(block_in_nodes, seqs_id)
                self.add_node(seqs_id=seqs_id,
                              id=current_node_id,
                              base=nucl,
                              in_nodes=in_nodes,
                              aligned_to=PangraphBuilderFromDAG.get_next_aligned_node_id(i, column_nodes_ids))
                PangraphBuilderFromDAG.update_block_in_nodes(block_in_nodes, seqs_id, current_node_id)

        self.add_block_out_edges_to_free_edges(block_in_nodes, block)
        self.manage_endings(block_in_nodes, block)

    @staticmethod
    def get_nodes_count(dagmaf: DAGMaf) -> int:
        nodes_count = 0
        SeqMafInfo = namedtuple('SeqMafInfo', ['maf_pos_count', 'srcSize'])
        seq_id_to_seqmafinfo = {}

        for n in dagmaf.dagmafnodes:
            # update seq_info
            for s in n.alignment:
                if s.id in seq_id_to_seqmafinfo:
                    s_info = seq_id_to_seqmafinfo[s.id]
                    seq_id_to_seqmafinfo[s.id] = SeqMafInfo(maf_pos_count=s.annotations["size"] + s_info.maf_pos_count,
                                                            srcSize=s_info.srcSize)
                else:
                    seq_id_to_seqmafinfo[s.id] = SeqMafInfo(maf_pos_count=s.annotations["size"],
                                                            srcSize=s.annotations["srcSize"])
            # sum number of current nodes
            number_of_columns = len(n.alignment[0].seq)
            for col_id in range(number_of_columns):
                letters_in_columns = set([n.alignment[i].seq[col_id] for i in range(len(n.alignment))]).difference(set('-'))
                nodes_count += len(letters_in_columns)

        nodes_to_complement_count = sum([info.srcSize - info.maf_pos_count for s, info in seq_id_to_seqmafinfo.items()])
        return nodes_count + nodes_to_complement_count

    @staticmethod
    def get_next_aligned_node_id(current_column_i, column_nodes_ids):
        if len(column_nodes_ids) > 1:
            return column_nodes_ids[(current_column_i + 1) % len(column_nodes_ids)]
        return None

    @staticmethod
    def get_starting_blocks(dagmaf):
        blocks_ids = [node.id for node in dagmaf.dagmafnodes]
        blocks_targeted_by_any_edge = [edge.to for node in dagmaf.dagmafnodes for edge in node.out_edges]
        starting_blocks_ids = set(blocks_ids)-set(blocks_targeted_by_any_edge)
        return list(starting_blocks_ids)



