from Bio import AlignIO
from itertools import islice
import numpy as np

import toolkit as toolkit
import nucleotides as nucleotides
from POAGraph import POAGraph
from Sequence import Source
from Node import Node


def parse_to_poagraphs(file_path, merge_option, multialignment_name, output_dir):
    maf_blocks = [*AlignIO.parse(file_path, "maf")]
    block_to_merge_ranges = _parse_merge_option_to_ranges(merge_option, len(maf_blocks))

    poagraphs = []
    for i, r in enumerate(block_to_merge_ranges):
        current_range_blocks = [*islice(maf_blocks, r.start, r.stop)]
        poagraph = POAGraph(name=multialignment_name,
                            title=multialignment_name + '_' + str(i),
                            path=toolkit.create_child_dir(output_dir, multialignment_name + '_' + str(i)),
                            version='NOVEMBER')
        poagraph = _blocks_to_poagraph(poagraph, current_range_blocks)
        poagraphs.append(poagraph)
    return poagraphs


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
                new_source = Source(ID=len(sources), name=line.name, title=line.name)
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


def get_blocks(file_path, multialignment_name, output_dir):
    #todo to powinno być jakoś zmergowane z wczytywaniem multialignmentu do poagrafów
    def find_next_block(src_name, next_start, maf_blocks):
        for i, block in enumerate(maf_blocks):
            for seqrec in block:
                if seqrec.name == src_name:
                    if seqrec.annotations["start"] == next_start:
                        return i
        return None

    maf_blocks = [*AlignIO.parse(file_path, "maf")]
    sources, sources_name_to_ID = _get_sources(maf_blocks)
    blocks = [Block() for i in range(len(maf_blocks))]
    for i, block in enumerate(maf_blocks):
        for seqrec in block:
            blocks[i].ID = i
            start = seqrec.annotations["start"]
            next_start = start + seqrec.annotations["size"]
            if next_start == seqrec.annotations["srcSize"]:
                blocks[i].srcID_to_next_blockID[sources_name_to_ID[seqrec.name]] = None
                continue
            next_block_ID = find_next_block(seqrec.name, next_start, maf_blocks)
            blocks[i].srcID_to_next_blockID[sources_name_to_ID[seqrec.name]] = next_block_ID
            blocks[i].srcID_to_strand[sources_name_to_ID[seqrec.name]] = seqrec.annotations["strand"]
    return blocks



class Block(object):
    def __init__(self):
        self.ID = -1
        self.srcID_to_next_blockID = {}
        self.srcID_to_strand = {}

    def __str__(self):
        return "ID: {0}\nsrcID_to_next_blockID: {1}\nsrcID_to_strand: {2}".format(self.ID, self.srcID_to_next_blockID, self.srcID_to_strand)
# for seqrec in multiple_alignment:
#     print("starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
#           (seqrec.annotations["start"],
#            seqrec.annotations["strand"],
#            seqrec.annotations["srcSize"],
#            seqrec.annotations["size"]))
