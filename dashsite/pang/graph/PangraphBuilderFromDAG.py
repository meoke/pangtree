import sys
from collections import deque, namedtuple

from fileformats.maf.DAGMaf import DAGMaf
from graph import nucleotides
from graph.FastaSource import FastaSource
from graph.Node import Node
from graph.PangraphBuilder import PangraphBuilder
from metadata.MultialignmentMetadata import MultialignmentMetadata


class EmptyIn:
    def __init__(self, to_block, node_id, seq_id, seq_pos):
        self.to_block = to_block
        self.node_id = node_id
        self.seq_id = seq_id
        self.seq_pos = seq_pos

    def __str__(self):
        return f"To: {self.to_block} " \
            f"Node_id: {self.node_id}, " \
            f"Seq_id: {self.seq_id}, Seq pos: {self.seq_pos}"


class EmptyOut:
    def __init__(self, from_block, node_id, to_block, seq_id, seq_pos, type):
        self.from_block = from_block
        self.node_id = node_id
        self.to_block = to_block
        self.seq_id = seq_id
        self.seq_pos = seq_pos
        self.type = type

    def __str__(self):
        return f"From: {self.from_block}, To: {self.to_block}, " \
            f"Type: {self.type}, Node_id: {self.node_id}, " \
            f"Seq_id: {self.seq_id}, Seq pos: {self.seq_pos}"


class SeqAttributes:
    def __init__(self, start, size, strand, srcSize):
        self.start = start
        self.size = size
        self.strand = strand
        self.srcSize = srcSize

    def get_last_pos(self):
        if self.strand == 1:
            return self.start + self.size - 1
        else:
            return self.srcSize - self.size + 1

    def __str__(self):
        return f"Start: {self.start}, Size: {self.size}"


class VisitOnlyOnceDeque():
    def __init__(self, iterable):
        self.d = deque(iterable)
        self.history = []

    def append(self, element):
        if element not in self.history:
            self.d.append(element)
            self.history.append(element)

    def popleft(self):
        element = self.d.popleft()
        self.history.append(element)
        return element

    def not_empty(self):
        if len(self.d):
            return True
        return False


class PangraphBuilderFromDAG(PangraphBuilder):
    def __init__(self, genomes_info: MultialignmentMetadata, fasta_source: FastaSource):
        super().__init__(genomes_info)
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

    def init_out_edges(self, dagmaf, pangraph):
        # todo it works for + strain only!
        PosInfo = namedtuple("PosInfo", ['start', 'block_id'])
        seq_id_to_first_pos = {seq_id: PosInfo(sys.maxsize, None) for seq_id in self.sequences_names}
        for n in dagmaf.dagmafnodes:
            for s in n.alignment:
                if s.annotations["start"] < seq_id_to_first_pos[s.id].start:
                    seq_id_to_first_pos[s.id] = PosInfo(s.annotations["start"], n.id)

        current_node_id = -1
        out_edges = []
        for seq_id, pos_info in seq_id_to_first_pos.items():
            in_node = []
            if pos_info.start == 0 or pos_info.start == sys.maxsize:
                continue
            else:
                for i in range(pos_info.start):
                    nucleotide = self.full_sequences[seq_id][i]
                    current_node_id += 1
                    pangraph._nodes[current_node_id] = Node(id=current_node_id,
                                                        base=nucleotides.code(nucleotide),
                                                        in_nodes=in_node,
                                                        aligned_to=None)
                    pangraph.add_path_to_node(path_name=seq_id, node_id=current_node_id)
                    in_node = [current_node_id]

                out_edges.append(EmptyOut(from_block=None,
                                          node_id=current_node_id,
                                          to_block=pos_info.block_id,
                                          seq_id=seq_id,
                                          seq_pos=pos_info.start-1,
                                          type=(1,-1),))
        return out_edges

    @staticmethod
    def get_max_node_id(pangraph):
        current_nodes_ids = [node.id for node in pangraph._nodes if node is not None]
        return max(current_nodes_ids) if current_nodes_ids else None

    def build(self, input, pangraph):
        sources_maf_info = self.get_sources_maf_info(input)
        nodes_count = PangraphBuilderFromDAG.get_nodes_count(input)
        pangraph._nodes = [None] * nodes_count
        pangraph._pathmanager.init_paths(self.sequences_names, nodes_count)


        out_edges = self.init_out_edges(input, pangraph)
        max_node_id = PangraphBuilderFromDAG.get_max_node_id(pangraph)
        current_node_id = max_node_id if max_node_id is not None else -1
        #to zbudować węzły, jeśli początku sekwencji nie ma w mafie (przejrzeć dag, znaleźć najmniejszy start i jeśłi jest większy niż 0 tzn że trzeba dociagac)
        #jesli byly budowane to dac id ostatniego do out edges, jesli nie to None
        in_edges = []
        blocks_deque = VisitOnlyOnceDeque(PangraphBuilderFromDAG.get_starting_blocks(input)) #todo czy pierwszy jest rzeczywiście pierwszy?
        while blocks_deque.not_empty():
            block_id = blocks_deque.popleft()
            block = input.dagmafnodes[block_id]
            block_width = len(block.alignment[0].seq)
            sequence_name_to_parameters = {seq.id: SeqAttributes(seq.annotations["start"], seq.annotations["size"], seq.annotations["strand"], seq.annotations["srcSize"]) for seq in block.alignment}
            local_edges = {seq.id : PangraphBuilderFromDAG.get_matching_connection(out_edges, seq.id, block_id, in_edges, sequence_name_to_parameters[seq.id]) for seq in block.alignment}
            #sprawdzac pozycje - jesli sie nie zgadzaja to robic dobudowki śródblokowe

            for col in range(block_width):
                sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block.alignment}
                nodes_codes = sorted([*(set([nucleotide for nucleotide in sequence_name_to_nucleotide.values()])).difference(set(['-']))])
                column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]

                for i, nucl in enumerate(nodes_codes):
                    current_node_id += 1
                    node = Node(id=current_node_id,
                                base=nucleotides.code(nucl),
                                in_nodes=set(),
                                aligned_to=PangraphBuilderFromDAG.get_next_aligned_node_id(i, column_nodes_ids))

                    for sequence, nucleotide in sequence_name_to_nucleotide.items():
                        if nucleotide == nucl:
                            pangraph.add_path_to_node(path_name=sequence, node_id=current_node_id)
                            # find previous node
                            last_node_id = local_edges[sequence]
                            if last_node_id is not None:
                                node.in_nodes.add(last_node_id)
                            elif sequence_name_to_parameters[sequence].start != 0:
                                in_edges.append(EmptyIn(block_id, current_node_id, sequence, sequence_name_to_parameters[sequence].start))
                            local_edges[sequence] = current_node_id
                    node.in_nodes = list(node.in_nodes)
                    pangraph._nodes[current_node_id] = node

            for edge in block.out_edges:
                if edge.edge_type == (1,-1):
                    blocks_deque.append(edge.to)
                for seq in edge.sequences:
                    seq_end = seq[1]
                    seq_id = seq_end.seq_id
                    out_edges.append(EmptyOut(block_id, local_edges[seq_id], edge.to, seq_id, sequence_name_to_parameters[seq_id].get_last_pos(), edge.edge_type))

            for seq in block.alignment:
                if sequence_name_to_parameters[seq.id].get_last_pos() != sequence_name_to_parameters[seq.id].srcSize:
                    add_out_edge = True
                for edge in block.out_edges:
                    edge_sequences_ids = [edge_seq[0].seq_id for edge_seq in edge.sequences]
                    if seq.id in edge_sequences_ids:
                        add_out_edge = False
                        break
                if add_out_edge:
                    out_edges.append(EmptyOut(block_id, local_edges[seq.id], None, seq.id, sequence_name_to_parameters[seq.id].get_last_pos(), (1,-1)))


        # do koncowek dociagnac wezly
        self.complement_endings(input, pangraph, in_edges)

        #todo krawędzie in-out możnaby od razu jakoś parować
        in_edge_to_remove = []
        out_edge_to_remove = []
        for i, out_edge in enumerate(out_edges):
            for j, in_edge in enumerate(in_edges):
                if out_edge.seq_pos == in_edge.seq_pos -1\
                        and out_edge.type == (1,-1): # todo it works for + strain only!
                    pangraph._nodes[in_edge.node_id].in_nodes.append(out_edge.node_id)
                    pangraph._nodes[in_edge.node_id].in_nodes = sorted(pangraph._nodes[in_edge.node_id].in_nodes)
                    in_edge_to_remove.append(j)
                    out_edge_to_remove.append(i)
                elif out_edge.seq_id == in_edge.seq_id and \
                        out_edge.to_block == in_edge.to_block and \
                        out_edge.type == (1, -1):
                    self.complement_nodes(pangraph, out_edge, in_edge)
                    out_edge_to_remove.append(i)
                    in_edge_to_remove.append(j)

        in_edges = [in_edge for i, in_edge in enumerate(in_edges) if i not in in_edge_to_remove]
        out_edges = [out_edge for i, out_edge in enumerate(out_edges) if i not in out_edge_to_remove]

        # for index in sorted(in_edge_to_remove, reverse=True):
        #     del in_edges[index]
        #
        # for index in sorted(out_edge_to_remove, reverse=True):
        #     del out_edges[index]
        # with open("a.txt", 'w') as o:
        #     for i in in_edges:
        #         o.write(f"{i}\n")
        #     for i in out_edges:
        #         o.write(f"{i}\n")

        pangraph._pathmanager.remove_empty_paths()

    def complement_endings(self, dagmaf, pangraph, in_edges):
        # todo it works for + strain only!
        #usuwac z out edges!
        PosInfo = namedtuple("PosInfo", ['start', 'size', 'srcSize', 'block_id'])
        seq_id_to_last_pos = {seq_id: PosInfo(-1, None, None, None) for seq_id in self.sequences_names}
        for n in dagmaf.dagmafnodes:
            for s in n.alignment:
                if s.annotations["start"] > seq_id_to_last_pos[s.id].start:
                    seq_id_to_last_pos[s.id] = PosInfo(s.annotations["start"], s.annotations["size"], s.annotations["srcSize"], n.id)

        current_node_id = PangraphBuilderFromDAG.get_max_node_id(pangraph) + 1
        for seq_id, pos_info in seq_id_to_last_pos.items():
            if pos_info.start == -1 or pos_info.start + pos_info.size == pos_info.srcSize:
                continue
            else:
                in_node=[]
                in_edges.append(EmptyIn(None, current_node_id, seq_id, pos_info.start+pos_info.size))
                for i in range(pos_info.start + pos_info.size, pos_info.srcSize):
                    nucleotide = self.full_sequences[seq_id][i]
                    pangraph._nodes[current_node_id] = Node(id=current_node_id,
                                                                base=nucleotides.code(nucleotide),
                                                                in_nodes=in_node,
                                                                aligned_to=None)
                    pangraph.add_path_to_node(path_name=seq_id, node_id=current_node_id)
                    in_node = [current_node_id]
                    current_node_id += 1

    def complement_nodes(self, pangraph, out_edge, in_edge):
        start_seq_pos = out_edge.seq_pos + 1
        end_seq_pos = in_edge.seq_pos
        current_node_id = PangraphBuilderFromDAG.get_max_node_id(pangraph) + 1
        seq_id = out_edge.seq_id
        in_node = [out_edge.node_id]
        for seq_pos in range(start_seq_pos, end_seq_pos):
            nucleotide = self.full_sequences[seq_id][seq_pos]
            pangraph._nodes[current_node_id] = Node(id=current_node_id,
                                                    base=nucleotides.code(nucleotide),
                                                    in_nodes=in_node,
                                                    aligned_to=None)
            pangraph.add_path_to_node(path_name=seq_id, node_id=current_node_id)
            in_node = [current_node_id]
            current_node_id += 1
        in_edge_in_nodes = pangraph._nodes[in_edge.node_id].in_nodes
        in_edge_in_nodes.append(current_node_id-1)
        pangraph._nodes[in_edge.node_id].in_nodes = sorted(in_edge_in_nodes)

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
            #sum number of current nodes
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
    def get_matching_connection(out_edges, seq_id, to_block_id, in_edges, sequence_params):
        out_edge_id = None
        in_edge_id = None
        edge = None
        for i, out_edge in enumerate(out_edges):
            if out_edge.seq_id == seq_id \
                    and out_edge.type == (1, -1) \
                    and out_edge.to_block == to_block_id \
                    and out_edge.seq_pos == sequence_params.start -1:
                out_edge_id = i
                edge = out_edge
                for j, in_edge in enumerate(in_edges):
                    if seq_id == in_edge.seq_id \
                            and to_block_id == in_edge.to_block \
                            and out_edge.seq_pos == in_edge.seq_pos - 1:
                        in_edge_id = j
        if in_edge_id is not None:
            del in_edges[in_edge_id]
        if out_edge_id is not None:
            del out_edges[out_edge_id]
            return edge.node_id
        return None

    def get_sources_maf_info(self, dagmaf):
        d = {}
        first_blocks = self.get_starting_blocks(dagmaf)
        dd = VisitOnlyOnceDeque

    @staticmethod
    def get_starting_blocks(dagmaf):
        blocks_ids = [node.id for node in dagmaf.dagmafnodes]
        blocks_targeted_by_any_edge = [edge.to for node in dagmaf.dagmafnodes for edge in node.out_edges]
        starting_blocks_ids = set(blocks_ids)-set(blocks_targeted_by_any_edge)
        return list(starting_blocks_ids)

