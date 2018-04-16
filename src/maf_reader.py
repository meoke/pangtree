from Bio import AlignIO
import numpy as np
from itertools import islice

import toolkit as t
import nucleotides as nucleotides
from Sequence import Source
from Node import Node
from copy import deepcopy
import networkx as nx
import random
from POAGraph import POAGraph
from Errors import NoSequenceContinuationFound


import matplotlib.pyplot as plt


def parse_to_poagraphs(file_path, merge_option, multialignment_name, output_dir):
    # statystyki

    # usunięcie cykli z bloków i zachowanie informacji, jak zbudować poagraf(y)
    mafblocks_nxgraph, sequences = _uncycle_blocks_version2(file_path)
    poagraphs = []

    # budowanie poagrafów
    for i, component_nodes in enumerate(nx.weakly_connected_components(mafblocks_nxgraph)):
        subgraph = mafblocks_nxgraph.subgraph(component_nodes)
        poagraph = POAGraph(name=multialignment_name,
                            title=multialignment_name + '_' + str(i),
                            path=t.create_child_dir(output_dir, multialignment_name + '_' + str(i)),
                            version='APRIL')
        poagraphs.append(_read_poagraph(poagraph, subgraph, sequences))

    return poagraphs
    #maf_blocks = [*AlignIO.parse(file_path, "maf")]

    #block_to_merge_ranges = _parse_merge_option_to_ranges(merge_option, len(maf_blocks))
    # poagraphs = []
    #todo wykorzystac wiedze z blocks do tworzenia poagraphs
    # for i, r in enumerate(block_to_merge_ranges):
    #     current_range_blocks = [*islice(maf_blocks, r.start, r.stop)]
    #     poagraph = POAGraph(name=multialignment_name + '_' + str(i),
    #                         title=multialignment_name + '_' + str(i),
    #                         path=t.create_child_dir(output_dir, multialignment_name + '_' + str(i)),
    #                         version='NOVEMBER')
    #     poagraph = _blocks_to_poagraph(poagraph, current_range_blocks)
    #     poagraph.set_sources_weights()
    #     poagraphs.append(poagraph)
    # return poagraphs, blocks


def _uncycle_blocks_version2(file_path: str) -> nx.DiGraph:
    def print_data(mafblocks_nxgraph):
        print("Edges count: {}".format(mafblocks_nxgraph.number_of_edges()))
        # print("Nodes count: {}".format(mafblocks_nxgraph.number_of_nodes()))

        # for node in mafblocks_nxgraph.nodes.data():
        #     print(node)
        # for edge in mafblocks_nxgraph.edges.data():
        #     print(edge)

        plt.subplot(121)
        nx.draw(mafblocks_nxgraph, with_labels=True, font_weight='bold')
        plt.show()

    def cycle_exist(cycles_iterator):
        try:
            first = next(cycles_iterator)
        except StopIteration:
            return False
        return True

    mafblocks_nxgraph = nx.DiGraph()
    sequences = {}

    # dodanie wszystkich bloków, bez krawędzi i bez sekwencji "Anc0"
    maf_blocks = AlignIO.parse(file_path, "maf")
    for i, block in enumerate(maf_blocks):
        seqrec_to_remove_ID = -1
        for j, seqrec in enumerate(block):
            short_seqname = seqrec.name.split(".")[0]
            if "Anc0" == short_seqname:
                seqrec_to_remove_ID = j
                continue
            try:
                sequences[short_seqname].append((i, seqrec.annotations["start"]))
            except KeyError:
                sequences[short_seqname] = [(i, seqrec.annotations["start"])]
        if seqrec_to_remove_ID != -1:
            del block._records[seqrec_to_remove_ID]
        mafblocks_nxgraph.add_node(i, mafblock=block)


    # sortowanie
    for seqname, blockid_start_list in sequences.items():
        sequences[seqname] = sorted(blockid_start_list, key=lambda blockid_start : blockid_start[1])

    # dodanie krawędzi, w zależności od przebiegu sekwencji (uproszczenie grafu)
    for node in mafblocks_nxgraph.nodes.data():
        for seqrec in node[1]["mafblock"]:
            short_seqname = seqrec.name.split(".")[0]

            try:
                next_blockID = _find_next_seqrec(node[0], seqrec, sequences)
                nextnode = mafblocks_nxgraph.nodes.data()[next_blockID]
                for seq in nextnode["mafblock"]:
                    if short_seqname in seq.name:
                        next_seqrec = seq
                        break
            except NoSequenceContinuationFound:
                continue

            if next_seqrec.annotations["strand"] == seqrec.annotations["strand"]:
                if (node[0], next_blockID) in mafblocks_nxgraph.edges:
                    mafblocks_nxgraph.edges[(node[0], next_blockID)]["weight"] += 1
                    mafblocks_nxgraph.edges[(node[0], next_blockID)]["seqs"].append(short_seqname)
                else:
                    mafblocks_nxgraph.add_edge(node[0], next_blockID, weight = 1, seqs=[short_seqname])

    print_data(mafblocks_nxgraph)
    # rozcięcie cykli (przejść krawędzie zaczynając od największej wagi i dodawać do finalnego grafu, jeśli
    # to nie spowoduje powstania cyklu)
    temp_dg = mafblocks_nxgraph.copy()
    all_edges = [(u, v) for (u, v) in temp_dg.edges()]
    temp_dg.remove_edges_from(all_edges)
    edges_causing_cycles = []
    for (u, v, wt) in sorted(mafblocks_nxgraph.edges.data('weight'), key=lambda x : x[2], reverse=True):
        temp_dg.add_edge(u, v, weight=wt)
        cycles = nx.simple_cycles(temp_dg)
        if cycle_exist(cycles):
            edges_causing_cycles.append((u, v))
            temp_dg.remove_edge(u, v)

    mafblocks_nxgraph.remove_edges_from(edges_causing_cycles)
    print_data(mafblocks_nxgraph)

    return mafblocks_nxgraph, sequences


def _find_next_seqrec(current_blockID: int, seqrec: AlignIO.MafIO.SeqRecord, sequences_sorted_by_pos: []) -> int:
    next_seqrec_start = seqrec.annotations["start"] + seqrec.annotations["size"]
    short_seqrecname = seqrec.name.split(".")[0]
    seqrec_info = sequences_sorted_by_pos[short_seqrecname]
    for i, blockID_start in enumerate(seqrec_info):
        if blockID_start[0] == current_blockID:
            try:
                next_blockID = seqrec_info[i+1][0]
            except Exception:
                raise NoSequenceContinuationFound("Sequence {} has no continuation after {} nucleotide.".format(
                    short_seqrecname, next_seqrec_start - 1))
    return next_blockID

    # for n in mafblocks_nxgraph.nodes.data():
    #     for s in n[1]["mafblock"]:
    #         short_sname = s.name.split(".")[0]
    #         if short_seqrecname == short_sname and s.annotations["start"] == next_seqrec_start:
    #             return n[0], s

#
# def _uncycle_blocks_version_1(file_path):
#     def cycle_exist(cycles_iterator):
#         try:
#             first = next(cycles_iterator)
#         except StopIteration:
#             return False
#         return True
#
#     def draw_blocks_graph(DG):
#         pos = nx.get_node_attributes(DG, 'pos')
#         nx.draw(DG, with_labels=True, font_weight='bold', pos=pos)
#         labels = nx.get_edge_attributes(DG, 'weight')
#         nx.draw_networkx_edge_labels(DG, pos=pos, edge_labels=labels)
#         plt.show()
#
#     def find_next_precise_block(src_name, next_start, maf_blocks):
#         for i, block in enumerate(maf_blocks):
#             for seqrec in block:
#                 if seqrec.name == src_name:
#                     if seqrec.annotations["start"] == next_start:
#                         return i
#         return None
#
#     def get_src_path(src_name, maf_blocks):
#         src_path = []
#         for i, block in enumerate(maf_blocks):
#             for seqrec in block:
#                 if seqrec.name == src_name:
#                     src_path.append((i, seqrec.annotations["start"]))
#                     break
#         return sorted(src_path, key=lambda blockID_start : blockID_start[1])
#
#     maf_blocks = [*AlignIO.parse(file_path, "maf")]
#
#     DG = nx.DiGraph()
#     DG.add_nodes_from(range(len(maf_blocks)))
#     for i, block in enumerate(maf_blocks):
#         rand_x = random.randrange(20)
#         rand_y = random.randrange(20)
#         DG.nodes[i]['pos'] = (rand_x, rand_y)
#
#     srcs, src_name_to_ID = _get_sources(maf_blocks)
#     for src in srcs:
#         src_path = get_src_path(src.name, maf_blocks)
#         for i, blockID_start in enumerate(src_path):
#             if i == len(src_path)-1:
#                 break
#             left_node = blockID_start[0]
#             right_node = (src_path[i+1])[0]
#             if right_node in list(DG.successors(left_node)):
#                 DG[left_node][right_node]['weight'] = DG[left_node][right_node]['weight'] + 1
#             else:
#                 DG.add_edge(left_node, right_node, weight=1, active=True)
#
#     temp_dg = DG.copy()
#     edges_to_remove = [(u, v) for (u, v) in temp_dg.edges()]
#     temp_dg.remove_edges_from(edges_to_remove)
#     for (u, v, wt) in sorted(DG.edges.data('weight'), key=lambda x : x[2], reverse=True):
#         temp_dg.add_edge(u, v, weight=wt)
#         cycles = nx.simple_cycles(temp_dg)
#         if cycle_exist(cycles):
#             DG.edges[u, v]['active'] = False
#             temp_dg.remove_edge(u, v)
#         else:
#             DG.edges[u, v]['active'] = True
#
#     return DG
#
#
# def _uncycle_blocks_simple_version(file_path, print_analysis=True):
#     def get_src_path(src_name, maf_blocks):
#         src_path = []
#         for i, block in enumerate(maf_blocks):
#             for seqrec in block:
#                 if seqrec.name == src_name:
#                     src_path.append((i, seqrec.annotations["start"]))
#                     break
#         return sorted(src_path, key=lambda blockID_start : blockID_start[1])
#
#     def print_blocks_analysis(maf_blocks):
#         blocks = {}
#         for i, block in enumerate(maf_blocks):
#             blocks[i] = [0, 0]
#             blocks_srcs = []
#             for line in block:
#                 blocks_srcs.append(line.name)
#                 if line.annotations["strand"] == 1:
#                     blocks[i][0] += 1
#                 elif line.annotations["strand"] == -1:
#                     blocks[i][1] += 1
#             minus = "MINUS" if blocks[i][1] > blocks[i][0] else ""
#             src = blocks_srcs[0] if len(blocks_srcs) == 1 else ""
#             print("block: {}, +: {}, -:{}, {}, {}".format(i, blocks[i][0], blocks[i][1], minus, src))
#
#     maf_blocks = [*AlignIO.parse(file_path, "maf")]
#     if print_analysis:
#         print_blocks_analysis(maf_blocks)
#
#
#     DG = nx.DiGraph()
#     DG.add_nodes_from(range(len(maf_blocks)))
#     srcs, src_name_to_ID = _get_sources(maf_blocks)
#
#     blocks = {}
#     for i, block in enumerate(maf_blocks):
#         DG.nodes[i]['pos'] = (random.randrange(20), random.randrange(20))
#
#         #trzeba określić czy jest większość +, jeśli nie, to obracamy
#         strands = [line.annotations["strand"] for line in block]
#         plus_strands_count = len([strand for strand in strands if strand == 1])
#         minus_strands_count = len([strand for strand in strands if strand == -1])
#
#         #obracanie
#         #uzupelnic
#
#         #dodanie tylko krawedzi +
#         # for line in block:
#         #     if line.annotations["strand"] == 1:
#         #         src = line.name
#         #         next_block_ID =
#
#         blocks[i] = [0,0]
#         blocks_srcs = []
#         for line in block:
#             blocks_srcs.append(line.name)
#             if line.annotations["strand"] == 1:
#                 blocks[i][0] += 1
#             elif line.annotations["strand"] == -1:
#                 blocks[i][1] += 1
#         minus = "MINUS" if blocks[i][1] > blocks[i][0] else ""
#         src = blocks_srcs[0] if len(blocks_srcs) == 1 else ""
#         print("block: {}, +: {}, -:{}, {}, {}".format(i, blocks[i][0], blocks[i][1], minus, src))
#
#     for src in srcs:
#         src_path = get_src_path(src.name, maf_blocks)
#         for i, blockID_start in enumerate(src_path):
#             if i == len(src_path)-1:
#                 break
#             left_node = blockID_start[0]
#             right_node = (src_path[i+1])[0]
#             if right_node in list(DG.successors(left_node)):
#                 DG[left_node][right_node]['weight'] = DG[left_node][right_node]['weight'] + 1
#             else:
#                 DG.add_edge(left_node, right_node, weight=1, active=True)
#

# class Block(object):
#     def __init__(self):
#         self.ID = -1
#         self.srcID_to_next_blockID = {}
#         self.srcID_to_strand = {}
#         self.next_blockID_to_weight = {}
#
#     def __str__(self):
#         return "ID: {0}\nsrcID_to_next_blockID: {1}\nsrcID_to_strand: {2}".format(self.ID, self.srcID_to_next_blockID, self.srcID_to_strand)

def _read_poagraph(poagraph, subgraph, sequences):
    def add_block(poagraph, nmatrix, piesek, currentnodeid, blockID, brakuje):
        column_nodes = {}
        source_to_its_last_node_ID = {source_ID: -1 for source_ID in range(nmatrix.shape[0])}
        sequences_nodes = {srcID: [] for srcID in range(nmatrix.shape[0])}
        nextnodeid = currentnodeid
        block_nodes = 0
        for column_ID in range(nmatrix.shape[1]):
            for row_ID in range(nmatrix.shape[0]):
                nucl = nmatrix[row_ID, column_ID]
                if nucl == nucleotides.code('-'):
                    continue

                if nucl not in column_nodes:
                    column_nodes[nucl]=Node(ID=nextnodeid,
                                            base=nucleotides.decode(nucl))
                    nextnodeid+=1
                    block_nodes +=1

                sequences_nodes[row_ID].append(column_nodes[nucl].ID)
                if source_to_its_last_node_ID[row_ID] != -1:
                    column_nodes[nucl].add_in_node([source_to_its_last_node_ID[row_ID]])
                source_to_its_last_node_ID[row_ID] = column_nodes[nucl].ID

                try:
                    poagraph.ns[row_ID, column_nodes[nucl].ID] = True
                except:
                    print("ech")


            sorted_column_nodes = sorted([*column_nodes.values()], key=lambda node: node.ID)
            for i, node in enumerate(sorted_column_nodes):
                if not len(sorted_column_nodes) == 1:
                    node.aligned_to = sorted_column_nodes[(i + 1) % len(column_nodes)].ID
                poagraph.nodes[node.ID] = node

            column_nodes = {}

        for srcid in piesek.keys():
            if sequences_nodes[srcid]:
                piesek[srcid][blockID] = (sequences_nodes[srcid][0], sequences_nodes[srcid][-1])
        return poagraph, block_nodes, piesek, brakuje

    # ?# z sequences odczytać, w którym bloku się kończy -> jeśli jest dziura, to zaciągać dane z fasty

    # współrzędne sekwencji:
    # jeśli w tym bloku się zaczyna - OK
    # jeśli wg sequences istenieje połączenie tego bloku i jakiegoś już dodanego bloku i w subgraphie takie połączenie też jest  - podpiąć
    # jeśli brakują nukleotydy, bo są nieujęte - zaciągnąć dane z fasty

    sources_name_to_ID = {}

    for id, s in enumerate(sorted(sequences.keys())):
        sources_name_to_ID[s] = id
        source = Source(ID=id, name=s, title=s)
        poagraph.add_source(source)

    # ustalić liczbę sekwencji
    all_sequences_count = len(sequences)

    # dla każdego bloku w subgraphie, obojętnie jaka kolejność
    # stworzyć macierz nukleotydów w nim
    blockid_to_nucleotides_matrix = {}
    for n in subgraph.nodes.data():
        blockid = n[0]
        mafblock = n[1]["mafblock"]
        # ustalić szerokość bloku (uwaga na same minusy w kolumnie?)
        block_width = len(mafblock[0].seq)

        # stworzyć macierz liczba sekwencji x szerokość bloku
        nmatrix = np.zeros(shape=(all_sequences_count, block_width), dtype=np.int8)

        # wypełnić macierz nukleotydami
        for i, line in enumerate(mafblock):
            # print("\r\t\tLine " + str(j+1) + '/' + str(len(block)), end='')
            short_sequence_name = line.name.split(".")[0]
            sequence_source_ID = sources_name_to_ID[short_sequence_name]
            for nucl_ID, nucl_base in enumerate(line.seq):
                nmatrix[sequence_source_ID][nucl_ID] = nucleotides.code(nucl_base)

        blockid_to_nucleotides_matrix[blockid] = nmatrix

    # sprawdzić ile będzie węzłów w sumie
    nodes_count = 0
    for nmatrix in blockid_to_nucleotides_matrix.values():
        for column_ID in range(nmatrix.shape[1]):
            nodes_count += len(set(nmatrix[:,column_ID]) - set([0]))

    # przygotować macierz połączeń węzłów i sekwencji
    print("All nodes count: {}".format(nodes_count))
    poagraph.ns = np.zeros(shape=(all_sequences_count, nodes_count), dtype=np.bool)
    print(nodes_count)

    # przygotować listę z miejscem na węzły
    poagraph.nodes = [None] * nodes_count

    # znalezc pierwsze bloki w subgraphie (bez ścieżek wchodzących, równoległe), żeby móc je obsłużyć jako pierwsze
    blocks_stack = []
    for blockid in subgraph.nodes:
        if subgraph.in_degree(blockid) == 0:
            blocks_stack.append(blockid)

    # każdy blok dodać do poagrafu przechodząc przez drzewo bloków depth-first-search
    # kolejnosć wyznaczają krawędzie obecne w
    processed_blocks = []
    piesek = {srcid:{} for srcid in range(len(sequences))}
    nextnodeid = 0
    brakuje = 0
    current_nodes_added_count = 0
    while True:
        next_blocks = []
        while blocks_stack:
            blockID_to_convert = blocks_stack.pop()
            processed_blocks.append(blockID_to_convert)
            print("procesuje {}".format(blockID_to_convert))
            poagraph, added_nodes, piesek, brakuje = add_block(poagraph, blockid_to_nucleotides_matrix[blockID_to_convert], piesek, current_nodes_added_count, blockID_to_convert, brakuje)
            current_nodes_added_count += added_nodes
            for start_end in subgraph.out_edges(blockID_to_convert):
                if start_end[1] not in processed_blocks and start_end[1] not in next_blocks and start_end[1] not in blocks_stack:
                    next_blocks.append(start_end[1])

        processed_blocks.extend(next_blocks)
        if next_blocks == []:
            break
        blocks_stack = next_blocks

    print(brakuje)
    # dodac krawedzie pomiedzy blokami, tam gdzie mozna (istnieją i są potrzebne)
    for srcname, srcblocks in sequences.items():
        srcid = sources_name_to_ID[srcname]
        for i in range(len(srcblocks)):
            if i < len(srcblocks)-1:
                if (srcblocks[i][0], srcblocks[i+1][0]) in subgraph.edges:
                    #id ostatniego wezla srcid w srcbklocks[i] i pierwszego w srcblocks[i+1]
                    right_node_id = piesek[srcid][srcblocks[i+1][0]][0]
                    left_node_id = piesek[srcid][srcblocks[i][0]][1]
                    try:
                        poagraph.nodes[right_node_id].add_in_node(left_node_id)
                    except:
                        pass

        print(s)

    poagraph.set_sources_weights()
    return poagraph
    # poagraph.set_sources_weights()

def _parse_merge_option_to_ranges(merge_option, blocks_count):
    if not merge_option:
        return [range(i, i+1) for i in range(blocks_count)]

    if merge_option == 'all':
        return [range(blocks_count)]

    else:
        def get_range(blocks_chunk):
            chunk_borders = blocks_chunk.split(':')
            if len(chunk_borders) == 1:
                return range(int(chunk_borders[0]), int(chunk_borders[0])+1)
            else:
                return range(int(chunk_borders[0]), int(chunk_borders[1])+1)

        return [get_range(blocks_chunk) for blocks_chunk in merge_option.split(',')]


def _blocks_to_poagraph(poagraph, blocks):
    all_blocks_width = _get_all_blocks_width(blocks)

    sources, sources_name_to_ID = _get_sources(blocks)
    for source in sources:
        poagraph.add_source(source)

    all_source_count = len(sources)
    nucleotides_matrix = np.zeros(shape=(all_blocks_width, all_source_count), dtype=np.int8) # zeros because '-' == 0

    # prepare nucleotides matrix
    current_column_ID = 0

    for i, block in enumerate(blocks):
        print("\tProcessing " + str(i+1) + "/" + str(len(blocks)) + " block.")  # todo logging
        for j, line in enumerate(block):
            print("\r\t\tLine " + str(j+1) + '/' + str(len(block)), end='')
            sequence_source_ID = sources_name_to_ID[line.name]  # todo logging
            for nucl_ID, nucl_base in enumerate(line.seq):
                nucleotides_matrix[current_column_ID + nucl_ID][sequence_source_ID] = nucleotides.code(nucl_base)

        current_column_ID += len(block._records[0])
        print('')

    #get nodes count
    nodes_count = 0
    for column_ID in range(nucleotides_matrix.shape[0]):
        nodes_count += len(set(nucleotides_matrix[column_ID][:]) - set([0]))

    #prepare matrix of connections between nodes and sources
    poagraph.ns = np.zeros(shape=(all_source_count, nodes_count), dtype=np.bool)

    #prepare place for nodes
    poagraph.nodes = [None] * nodes_count

    column_nodes = {}
    nodes_count = 0
    source_to_its_last_node_ID = {source_ID: -1 for source_ID, _ in enumerate(sources)}

    for column_ID in range(nucleotides_matrix.shape[0]):
        for row_ID in range(nucleotides_matrix.shape[1]):

            nucl = nucleotides_matrix[column_ID,row_ID]
            if nucl == nucleotides.code('-'):
                 continue

            if nucl not in column_nodes:
                column_nodes[nucl] = Node(ID=nodes_count, base=nucleotides.decode(nucl))
                nodes_count += 1

            node_id = column_nodes[nucl].ID

            if source_to_its_last_node_ID[row_ID] != -1:
                column_nodes[nucl].add_in_node([source_to_its_last_node_ID[row_ID]])
            source_to_its_last_node_ID[row_ID] = node_id

            try:
                poagraph.ns[row_ID, node_id] = True
            except:
                # pass
                raise ValueError("Wyjątek w maf readerze - kiedy to leci???")

            #progress = round(100 * (column_ID*nucleotides_matrix.shape[1]+row_ID+1) / all_nucles_count, 2) # todo logging
            #progress = column_ID #/ nucleotides_matrix.shape[0];
            #print("\r\tMultialignment ready in {0}".format(str(progress)), end='') # todo logging

        sorted_column_nodes = sorted([*column_nodes.values()], key=lambda node: node.ID)
        for i, node in enumerate(sorted_column_nodes):
            if not len(sorted_column_nodes) == 1:
                node.aligned_to = sorted_column_nodes[(i+1) % len(column_nodes)].ID
            poagraph.nodes[node.ID] = node

        column_nodes = {}
    print('') #todo logging
    return poagraph


def _get_sources(blocks):
    sources = []
    sources_name_to_ID = {}
    for block in blocks:
        for line in block:
            if line.name not in sources_name_to_ID:
                new_source = Source(ID=len(sources),
                                    name=line.name,
                                    title=line.name)
                sources.append(new_source)
                sources_name_to_ID[new_source.name] = new_source.ID
    return sources, sources_name_to_ID


def _get_all_blocks_width(blocks):
    return sum([len(block._records[0].seq) for block in blocks])


def _pretty_numpy_nucleotides_matrix_printer(numpy_matrix):
    for i in range(numpy_matrix.shape[1]):
        row = []
        for j in range(numpy_matrix.shape[0]):
            row.append(nucleotides.decode(numpy_matrix[j][i]))
        print("".join(row))


def _pretty_numpy_nucleotides_matrix_printer_normal(numpy_matrix):
    for i in range(numpy_matrix.shape[0]):
        line = []
        for j in range(numpy_matrix.shape[1]):
            line.append(nucleotides.decode(numpy_matrix[i][j]))
        print("".join(line))






# for seqrec in multiple_alignment:
#     print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
#           (seqrec.annotations["start"],
#            seqrec.annotations["strand"],
#            seqrec.annotations["srcSize"],
#            seqrec.annotations["size"]))
