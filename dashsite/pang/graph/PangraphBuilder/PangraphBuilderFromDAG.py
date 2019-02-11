import sys
from collections import deque, namedtuple

from mafgraph.graph import Block

from fileformats.maf.DAGMaf import DAGMaf
from graph import nucleotides
from graph.FastaSource import FastaSource
from graph.Node import Node
from graph.PangraphBuilder.PangraphBuilderBase import PangraphBuilderBase
from metadata.MultialignmentMetadata import MultialignmentMetadata
from mafgraph.mafreader import start_position

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


class PangraphBuilderFromDAG(PangraphBuilderBase):
    def __init__(self, genomes_info: MultialignmentMetadata, fasta_source: FastaSource):
        super().__init__(genomes_info)
        self.pangraph = None # będzie uzupełniany o węzły i połączenia
        self.dagmaf = None # przerabiany na pangraph
        # self.open_in = None # gromadzenie
        # self.open_out = None
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

    # def init_out_edges(self, dagmaf, pangraph):
    #     # todo it works for + strain only!
    #     PosInfo = namedtuple("PosInfo", ['start', 'block_id'])
    #     seq_id_to_first_pos = {seq_id: PosInfo(sys.maxsize, None) for seq_id in self.sequences_names}
    #     for n in dagmaf.dagmafnodes:
    #         for s in n.alignment:
    #             if s.annotations["start"] < seq_id_to_first_pos[s.id].start:
    #                 seq_id_to_first_pos[s.id] = PosInfo(s.annotations["start"], n.id)
    #
    #     current_node_id = -1
    #     out_edges = []
    #     for seq_id, pos_info in seq_id_to_first_pos.items():
    #         in_node = []
    #         if pos_info.start == 0 or pos_info.start == sys.maxsize:
    #             continue
    #         else:
    #             for i in range(pos_info.start):
    #                 nucleotide = self.full_sequences[seq_id][i]
    #                 current_node_id += 1
    #                 pangraph._nodes[current_node_id] = Node(id=current_node_id,
    #                                                     base=nucleotides.code(nucleotide),
    #                                                     in_nodes=in_node,
    #                                                     aligned_to=None)
    #                 pangraph.add_path_to_node(path_name=seq_id, node_id=current_node_id)
    #                 in_node = [current_node_id]
    #
    #             out_edges.append(EmptyOut(from_block=None,
    #                                       node_id=current_node_id,
    #                                       to_block=pos_info.block_id,
    #                                       seq_id=seq_id,
    #                                       seq_pos=pos_info.start-1,
    #                                       type=(1,-1),))
    #     return out_edges

    def get_max_node_id(self):
        current_nodes_ids = [node.id for node in self._nodes if node is not None]
        return max(current_nodes_ids) if current_nodes_ids else None

    def init_pangraph(self):
        nodes_count = PangraphBuilderFromDAG.get_nodes_count(self.input)
        self.pangraph._nodes = [None] * nodes_count
        self.pangraph._pathmanager.init_paths(self.sequences_names, nodes_count)

    # class BlockSequenceRecord:
    #     def __init__(self, start_node_id, end_node_id):
    #         self.start_node_id = start_node_id
    #         self.end_node_id = end_node_id

    SInfo = namedtuple('SInfo', ['block_id',
                                 'start',
                                 'strand',
                                 'size',
                                 'srcSize'])

    Edge = namedtuple('Edge', ['seq_id', 'from_block_id', 'last_node_id', 'to_block_id'])

    def build(self, input, pangraph):
        self.dagmaf = input
        self.pangraph = pangraph
        self.init_pangraph()

        self.free_edges = {seq_id: [] for seq_id in self.sequences_names}
        # self.blocks_info = {block.id: dict() for block in self.dagmaf.dagmafnodes}
        self.set_seqs_info()

        self.complement_starting_nodes()
        for mafnode in self.dagmaf.dagmafnodes:
            self.process_block(mafnode)
        self.complement_ending_nodes()

    def set_seqs_info(self):
        self.seqs_info = {seq_id: [] for seq_id in self.sequences_names}

        for n in self.dagmaf.dagmafnodes:
            for seq in n:
                self.seqs_info[seq.id].append(PangraphBuilderFromDAG.SInfo(block_id=n.id,
                                                                           start=start_position(seq),
                                                                           strand=seq.annotations["strand"],
                                                                           size=seq.annotations["size"],
                                                                           srcSize=seq.annotations["srcSize"]))
        for seq_id, seq_info_list in self.seqs_info.items():
            self.seqs_info[seq_id] = sorted(seq_info_list, key=lambda si: si.start)

    def get_left_sinfo(self, edge, seq_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == edge.from_block_id:
                return sinfo

    def get_edge_sinfos(self, edge,seq_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == edge.from_block_id:
                left_seq_info = sinfo
            if sinfo.block_id == edge.to_block_id:
                right_seq_info = sinfo
        return left_seq_info, right_seq_info

    def complement_starting_nodes(self):
        for seq_id, seq_info_list in self.seqs_info.items():
            first_pos = seq_info_list[0].start
            if first_pos != 0:
                self.complement_sequence_starting_nodes(seq_id, first_pos)

    def complement_ending_nodes(self):
        for seq_id, edges in self.free_edges.items():
            if len(edges) == 1:
                last_edge = edges[0]
                self.complement_sequence_ending_nodes(seq_id, last_edge)
            if len(edges) > 1:
                raise Exception("Why there are many free edges?")

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
                                                                   from_block=None,
                                                                   last_node_id=current_node_id,
                                                                   to_block=first_block_sinfo.block_id))

    def complement_sequence_ending_nodes(self, seq_id, last_edge: Edge):
        current_node_id = self.get_max_node_id()
        in_nodes = [last_edge.last_node_id]
        last_block_sinfo = self.get_sinfo(last_edge)
        for i in range(last_block_sinfo.start + last_block_sinfo.size, 0):
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
                                        base=base,
                                        in_nodes=in_nodes,
                                        aligned_to=aligned_to)
        for seq_id in seqs_id:
            self.pangraph.add_path_to_node(path_name=seq_id, node_id=id)

    def get_block_in_nodes(self, block):
        edges_ids_to_remove = []
        block_in_nodes = dict()
        for seq in block.alignment:
            for i, edge in enumerate(self.free_edges[seq.id]):
                if edge.to_block_id == block.id:
                    block_in_nodes[seq.id] = edge.last_node_id
                    edges_ids_to_remove.append(i)
                    continue
            block_in_nodes[seq.id] = None
        self.free_edges = [edge for i, edge in enumerate(self.free_edges) if i not in edges_ids_to_remove]
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
        return [node_id for seq_id, node_id in block_in_nodes.items() if seq_id in seqs_id]

    @staticmethod
    def update_block_in_nodes(block_in_nodes, seqs_id, current_node_id):
        for seq_id in seqs_id:
            block_in_nodes[seq_id] = current_node_id

    def complement_sequence_middles(self, block_in_nodes):
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
                                                                   from_block=None,
                                                                   last_node_id=current_node_id,
                                                                   to_block=first_block_sinfo.block_id))

    def get_sinfo(self, seq_id, block_id):
        for sinfo in self.seqs_info[seq_id]:
            if sinfo.block_id == block_id:
                return sinfo

    @staticmethod
    def continuous_sequence(left_block_sinfo, right_block_sinfo):
        return left_block_sinfo.start + left_block_sinfo.size == right_block_sinfo.start

    def complement_sequence_middle_nodes(self, seq_id, last_pos, next_pos) -> int:
        return 0

    def add_block_out_edges_to_free_edges(self, block_in_nodes, block: Block):
        extended_sequences = []
        for edge in block.out_edges:
            for seq in edge.sequences:
                extended_sequences.append(seq[0].id)  # nie jestem pewna czy to dobrze (jak to jest w mafgraphie?)
            if edge.type != (1, -1):
                continue
            for seq in edge.sequences:
                left_block_sinfo, right_block_sinfo = self.get_edge_sinfos(edge, seq.id)
                if not PangraphBuilderFromDAG.continuous_sequence(left_block_sinfo, right_block_sinfo):
                    last_node_id = self.complement_sequence_middle_nodes(seq_id=seq.id,
                                                          last_pos=left_block_sinfo.start + left_block_sinfo.size,
                                                          next_pos = right_block_sinfo.start)
                    from_block_id = None
                else:
                    last_node_id = block_in_nodes[seq.id]
                    from_block_id = block.id

                self.free_edges[seq].append(
                    PangraphBuilderFromDAG.Edge(
                        seq_id=seq,
                        from_block_id=from_block_id,
                        last_node_id=last_node_id,
                        to_block_id=edge.to))
        for seq in block.alignment:
            last_block_sinfo: PangraphBuilderFromDAG.SInfo = self.get_sinfo(seq.id, block.id)
            if seq.id not in extended_sequences:
                if last_block_sinfo.start + last_block_sinfo.size != last_block_sinfo.srcSize:
                    last_node_id = self.complement_sequence_middle_nodes(seq_id=seq.id,
                                                                         last_pos=last_block_sinfo.start + last_block_sinfo.size,
                                                                         next_pos=last_block_sinfo.srcSize)
                    from_block_id = None
                else:
                    last_node_id = block_in_nodes[seq.id]
                    from_block_id = block.id
                self.free_edges[seq.id].append(PangraphBuilderFromDAG.Edge(seq_id=seq.id,
                                                                           from_block_id=from_block_id,
                                                                           last_node_id=last_node_id,
                                                                           to_block_id=None))

    def process_block(self, block):
        current_node_id = self.get_max_node_id()
        block_width = len(block.alignment[0].seq)
        sequence_name_to_parameters = {
        seq.id: PangraphBuilderFromDAG.SInfo(block_id=block.id,
                                             start=seq.annotations["start"],
                                             strand=seq.annotations["strand"],
                                             size=seq.annotations["size"],
                                             srcSize=seq.annotations["srcSize"]) for seq in block.alignment}

        block_in_nodes = self.get_block_in_nodes(block)

        for col in range(block_width):
            sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block.alignment}
            nodes_codes = PangraphBuilderFromDAG.get_column_nucleotides_codes(sequence_name_to_nucleotide)
            column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]
            for i, nucl in enumerate(nodes_codes):
                current_node_id += 1
                seqs_id = [seq_id for seq_id, n in sequence_name_to_nucleotide.items() if n == nucl]
                in_nodes = self.get_in_nodes()
                self.add_node(seqs_id=seqs_id,
                              id=current_node_id,
                              base=nucleotides.code(nucl),
                              in_nodes=in_nodes,
                              aligned_to=PangraphBuilderFromDAG.get_next_aligned_node_id(i, column_nodes_ids))
                PangraphBuilderFromDAG.update_block_in_nodes(block_in_nodes, seqs_id, current_node_id)

        self.add_block_out_edges_to_free_edges(block_in_nodes, block)




    # przy procesowaniu bloku
    # sprawdzić w seqs_info czy dany blok jest jako ostatnio w seqs_info dla tej sekwencji
    # jeśli nie - to nic
    # jeśli tak - dorzucić krawędź z tego bloku do None do edges

    # for seq_id, seq_info_list in self.seqs_info.items():
    #     last_pos = seq_info_list[-1].start + seq_info_list[-1].size
    #     if last_pos != seq_info_list[0].srcSize and self.exists_edge(seq_id):
    #         self.complement_sequence_ending_nodes(seq_id, last_pos)
    #         self.complement_sequence_ending_nodes(seq_id, last_pos)


    def build2(self, input, pangraph):
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

    # def complement_endings(self, dagmaf, pangraph, in_edges):
    #     # todo it works for + strain only!
    #     #usuwac z out edges!
    #     PosInfo = namedtuple("PosInfo", ['start', 'size', 'srcSize', 'block_id'])
    #     seq_id_to_last_pos = {seq_id: PosInfo(-1, None, None, None) for seq_id in self.sequences_names}
    #     for n in dagmaf.dagmafnodes:
    #         for s in n.alignment:
    #             if s.annotations["start"] > seq_id_to_last_pos[s.id].start:
    #                 seq_id_to_last_pos[s.id] = PosInfo(s.annotations["start"], s.annotations["size"], s.annotations["srcSize"], n.id)
    #
    #     current_node_id = PangraphBuilderFromDAG.get_max_node_id(pangraph) + 1
    #     for seq_id, pos_info in seq_id_to_last_pos.items():
    #         if pos_info.start == -1 or pos_info.start + pos_info.size == pos_info.srcSize:
    #             continue
    #         else:
    #             in_node=[]
    #             in_edges.append(EmptyIn(None, current_node_id, seq_id, pos_info.start+pos_info.size))
    #             for i in range(pos_info.start + pos_info.size, pos_info.srcSize):
    #                 nucleotide = self.full_sequences[seq_id][i]
    #                 pangraph._nodes[current_node_id] = Node(id=current_node_id,
    #                                                             base=nucleotides.code(nucleotide),
    #                                                             in_nodes=in_node,
    #                                                             aligned_to=None)
    #                 pangraph.add_path_to_node(path_name=seq_id, node_id=current_node_id)
    #                 in_node = [current_node_id]
    #                 current_node_id += 1
    #
    # def complement_nodes(self, pangraph, out_edge, in_edge):
    #     start_seq_pos = out_edge.seq_pos + 1
    #     end_seq_pos = in_edge.seq_pos
    #     current_node_id = PangraphBuilderFromDAG.get_max_node_id(pangraph) + 1
    #     seq_id = out_edge.seq_id
    #     in_node = [out_edge.node_id]
    #     for seq_pos in range(start_seq_pos, end_seq_pos):
    #         nucleotide = self.full_sequences[seq_id][seq_pos]
    #         pangraph._nodes[current_node_id] = Node(id=current_node_id,
    #                                                 base=nucleotides.code(nucleotide),
    #                                                 in_nodes=in_node,
    #                                                 aligned_to=None)
    #         pangraph.add_path_to_node(path_name=seq_id, node_id=current_node_id)
    #         in_node = [current_node_id]
    #         current_node_id += 1
    #     in_edge_in_nodes = pangraph._nodes[in_edge.node_id].in_nodes
    #     in_edge_in_nodes.append(current_node_id-1)
    #     pangraph._nodes[in_edge.node_id].in_nodes = sorted(in_edge_in_nodes)

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

    # @staticmethod
    # def get_matching_connection(out_edges, seq_id, to_block_id, in_edges, sequence_params):
    #     out_edge_id = None
    #     in_edge_id = None
    #     edge = None
    #     for i, out_edge in enumerate(out_edges):
    #         if out_edge.seq_id == seq_id \
    #                 and out_edge.type == (1, -1) \
    #                 and out_edge.to_block == to_block_id \
    #                 and out_edge.seq_pos == sequence_params.start -1:
    #             out_edge_id = i
    #             edge = out_edge
    #             for j, in_edge in enumerate(in_edges):
    #                 if seq_id == in_edge.seq_id \
    #                         and to_block_id == in_edge.to_block \
    #                         and out_edge.seq_pos == in_edge.seq_pos - 1:
    #                     in_edge_id = j
    #     if in_edge_id is not None:
    #         del in_edges[in_edge_id]
    #     if out_edge_id is not None:
    #         del out_edges[out_edge_id]
    #         return edge.node_id
    #     return None
    #
    # def get_sources_maf_info(self, dagmaf):
    #     d = {}
    #     first_blocks = self.get_starting_blocks(dagmaf)
    #     dd = VisitOnlyOnceDeque

    @staticmethod
    def get_starting_blocks(dagmaf):
        blocks_ids = [node.id for node in dagmaf.dagmafnodes]
        blocks_targeted_by_any_edge = [edge.to for node in dagmaf.dagmafnodes for edge in node.out_edges]
        starting_blocks_ids = set(blocks_ids)-set(blocks_targeted_by_any_edge)
        return list(starting_blocks_ids)



