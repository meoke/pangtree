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

    sources, sources_name_to_ID = _get_sources(blocks, all_blocks_width)
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

    column_nodes = {}
    nodes_count = 0
    source_to_its_last_node_ID = {source_ID: -1 for source_ID, _ in enumerate(sources)}
    all_nucles_count = nucleotides_matrix.shape[0]*nucleotides_matrix.shape[1]
    print(nucleotides_matrix.shape[0])
    print('')
    for column_ID in range(nucleotides_matrix.shape[0]):
        for row_ID in range(nucleotides_matrix.shape[1]):

            nucl = nucleotides_matrix[column_ID,row_ID]
            if nucl == nucleotides.code('-'):
                 continue

            if nucl not in column_nodes:
                column_nodes[nucl] = Node(ID=nodes_count, base=nucleotides.decode(nucl))  #todo usunac sources?
                nodes_count += 1

            # todo gdyby usunac sources, to ponizsze mozna zastapic powyzszym
            # if nucl not in column_nodes:
            #      column_nodes[nucl] = Node(ID=nodes_count, base=nucleotides.decode(nucl), sources=set([row_ID]))  #todo usunac sources?
            #      nodes_count += 1
            # else:
            #     column_nodes[nucl].sources.add(row_ID)
            node_id = column_nodes[nucl].ID

            if source_to_its_last_node_ID[row_ID] != -1:
                column_nodes[nucl].add_in_node([source_to_its_last_node_ID[row_ID]])
                #column_nodes[nucl].in_nodes.update([source_to_its_last_node_ID[row_ID]])
            source_to_its_last_node_ID[row_ID] = node_id
            #_update_source_sequence_info(row_ID, node_id, sources)
            try:
                poagraph.sources[row_ID].nodes_IDs[sources[row_ID].nodes_count] = node_id
                poagraph.sources[row_ID].nodes_count += 1
            except:
                pass

            # todo bo teraz aligned nodes to nie set
            # for key, node in column_nodes.items():
            #     other_nodes = column_nodes.copy()
            #     del other_nodes[key]
            #     node.aligned_to.update([node.currentID for node in other_nodes.values()])
            #     poagraph.add_node(node)

            #progress = round(100 * (column_ID*nucleotides_matrix.shape[1]+row_ID+1) / all_nucles_count, 2) # todo logging
            #progress = column_ID #/ nucleotides_matrix.shape[0];
            #print("\r\tMultialignment ready in {0}".format(str(progress)), end='') # todo logging

        sorted_column_nodes = sorted([*column_nodes.values()], key=lambda node: node.ID)
        for i, node in enumerate(sorted_column_nodes):
            if not len(sorted_column_nodes) == 1:
                node.aligned_to = sorted_column_nodes[(i+1) % len(column_nodes)].ID
            poagraph.add_node(node)

        column_nodes = {}
    poagraph.clean()
    print('') #todo logging
    return poagraph


def _get_sources(blocks, max_nodes_count):
    sources = []
    sources_name_to_ID = {}
    for block in blocks:
        for line in block:
            if line.name not in sources_name_to_ID:
                new_source = Source(ID=len(sources), name=line.name, title=line.name, max_nodes_count=max_nodes_count)
                sources.append(new_source)
                sources_name_to_ID[new_source.name] = new_source.ID
    return sources, sources_name_to_ID
#
#
def _get_all_blocks_width(blocks):
    return sum([len(block._records[0].seq) for block in blocks])# for lambda block : len(block._records[0].seq))
    # width = 0
    # for block in blocks:
    #     width += len(block._records[0].seq)  # all records in one block have the same length because they are aligned
    # return width
#
#
# def _update_source_sequence_info(source_ID, node_ID, sources):
#     sources[source_ID].add_node_ID(node_ID)
#
#
# def _pretty_numpy_nucleotides_matrix_printer(numpy_matrix):
#     for i in range(numpy_matrix.shape[1]):
#         row = []
#         for j in range(numpy_matrix.shape[0]):
#             row.append(nucleotides.decode(numpy_matrix[j][i]))
#         print("".join(row))