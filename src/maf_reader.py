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
    block_to_merge_ranges = _prepare_merge_ranges(merge_option, len(maf_blocks))

    poagraphs  = []
    for i, r in enumerate(block_to_merge_ranges):
        current_range_blocks = [*islice(maf_blocks, r.start, r.stop)]
        poagraph = POAGraph(name = multialignment_name,
                            title = multialignment_name + '_' + str(i),
                            path = toolkit.create_child_dir(output_dir, multialignment_name + '_' + str(i)),
                            version = 'NOVEMBER')
        poagraph = _process_blocks_for_poagraph(poagraph, current_range_blocks)
        poagraphs.append(poagraph)
    return poagraphs


def _prepare_merge_ranges(merge_option, blocks_count):
    if not merge_option:
        return [range(i, i+1) for i in range(blocks_count)]

    if merge_option == 'all':
        return [range(0, blocks_count)]

    else:
        def get_range(blocks_chunk):
            chunk_borders = blocks_chunk.split(':')
            if len(chunk_borders) == 1:
                return range(int(chunk_borders[0]), int(chunk_borders[0])+1)
            else:
                return range(int(chunk_borders[0]), int(chunk_borders[1])+1)

        return [get_range(blocks_chunk) for blocks_chunk in merge_option.split(',')]


def _process_blocks_for_poagraph(poagraph, blocks):
    sources, sources_name_to_ID = _extract_sources_info(blocks)

    for source in sources:
        poagraph.add_source(source)

    all_blocks_width = _get_all_blocks_width(blocks)
    all_source_count = len(sources)
    sequences = np.zeros(shape=(all_blocks_width, all_source_count), dtype=np.int8) # zeros because '-' == 0

    current_column_ID = 0
    for block in blocks:
        for line in block:
            sequence_source_ID = sources_name_to_ID[line.name]
            for nucl_ID, nucl_base in enumerate(line.seq):
                sequences[current_column_ID + nucl_ID][sequence_source_ID] = nucleotides.code(nucl_base)
        current_column_ID += len(block._records[0])

    column_nodes = {}
    nodes_count = 0
    source_to_its_last_node_ID = {source_ID: -1 for source_ID, _ in enumerate(sources)}
    for column_ID in range(sequences.shape[0]):
        for row_ID in range(sequences.shape[1]):

            nucl = sequences[column_ID,row_ID]
            if nucl == nucleotides.code('-'):
                 continue

            if nucl not in column_nodes:
                 column_nodes[nucl] = Node(ID = nodes_count, base=nucleotides.decode(nucl))
                 nodes_count += 1
            node_id = column_nodes[nucl].ID
            column_nodes[nucl].sources_count += 1

            if source_to_its_last_node_ID[row_ID] != -1:
                 column_nodes[nucl].in_nodes.update([source_to_its_last_node_ID[row_ID]])
            source_to_its_last_node_ID[row_ID] = node_id
            _update_source_sequence_info(row_ID, node_id, sources)

            for key, node in column_nodes.items():
                other_nodes = column_nodes.copy()
                del other_nodes[key]
                node.aligned_to.update([node.ID for node in other_nodes.values()])
                poagraph.add_node(node)
        column_nodes = {}

    return(poagraph)


def _extract_sources_info(blocks):
    sources = []
    sources_name_to_ID = {}
    for block in blocks:
        for line in block:
            if line.name not in sources_name_to_ID:
                new_source = Source(ID = len(sources), name = line.name, title = line.name)
                sources.append(new_source)
                sources_name_to_ID[new_source.name] = new_source.ID
    return (sources, sources_name_to_ID)


def _get_all_blocks_width(blocks):
    width = 0
    for block in blocks:
        width += len(block._records[0].seq) #all records have the same length
    return width


def _update_source_sequence_info(source_ID, node_ID, sources):
    if sources[source_ID].first_node_id == -1:
        sources[source_ID].first_node_id = node_ID
    #sources[source_ID].add_node_ID(node_ID)

    (sources[source_ID]).nodes_IDs.add(node_ID)
    pass


def _pretty_numpy_nucletides_matrix_printer(self, numpy_matrix):
    for i in range(numpy_matrix.shape[1]):
        row = []
        for j in range(numpy_matrix.shape[0]):
            row.append(nucleotides.decode(numpy_matrix[j][i]))
        print("".join(row))