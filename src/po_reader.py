import re
import numpy as np
from POAGraph import POAGraph
from Sequence import Consensus, Source
from Node import Node
# from NoConsensusFound import *


def parse_to_poagraph(file_path, output_dir):
    print('\tBuliding poagraph from ' + file_path) # todo logging
    with open(file_path) as po:
        #po_lines = po.readlines()
        poagraph = POAGraph( version=_read_value(po.readline()),#po_lines[0])
                             name=_read_value(po.readline()),#po_lines[1]),
                             title=_read_value(po.readline()),#po_lines[2]),
                             path=output_dir)
        poagraph.sources, poagraph.consensuses = _read_sequence_info(po)
        fill_nodes_info(poagraph, po) #uwaga, czy po jest zupatedowane?
        #_read_nodes_from_po_lines(poagraph, po_lines, len(poagraph.sources)-1)
    return poagraph


def _read_value(line):
    return line.split('=')[1].strip()


def _read_sequence_info(po_file_handler):
    length = int(_read_value(po_file_handler.readline()))
    source_count = int(_read_value(po_file_handler.readline()))

    source_ID = -1
    consensus_ID = -1
    sources = []
    consensuses = []
    for line in po_file_handler:

        sequence_name = _read_value(line)
        detailed_info_line = po_file_handler.readline()
        detailed_info = _read_value(detailed_info_line).split(' ')
        if 'CONSENS' in sequence_name:
            consensus_ID += 1
            consensus = Consensus(ID=consensus_ID,
                                  name=sequence_name,
                                  title=" ".join(detailed_info[4:]),
                                  max_nodes_count=int(detailed_info[0]))
            consensuses.append(consensus)
        else:
            source_ID += 1
            source = Source(ID=source_ID,
                            name=sequence_name,
                            title=" ".join(detailed_info[4:]),
                            weight=int(detailed_info[2]),
                            max_nodes_count=int(detailed_info[0]))
            #source.consensusID = int(source_info[3])#todo usuniete consensus ID, wiec tracimy ta informacje...
            sources.append(source)

        if source_ID + consensus_ID + 2 == source_count:
            break

    return sources, consensuses

def _read_node_parameters(node, code_letter):
    pattern = '{0}\d+'.format(code_letter)
    values_with_prefix_letters = re.findall(pattern, node)
    return [int(letter_value[1:]) for letter_value in values_with_prefix_letters]

def fill_nodes_info(poagraph, po_file_handler):
    def assign_this_node_to_its_sequences(sequences_IDs, node_ID):
        poagraph_sources_count = len(poagraph.sources)
        for sequence_ID in sequences_IDs:
            if sequence_ID < poagraph_sources_count:
                source_nodes_count = poagraph.sources[sequence_ID].nodes_count
                poagraph.sources[sequence_ID].nodes_IDs[source_nodes_count] = node_ID
                poagraph.sources[sequence_ID].nodes_count += 1
            else:
                consensus_ID = sequence_ID - poagraph_sources_count
                consensus_nodes_count = poagraph.consensuses[consensus_ID].nodes_count
                poagraph.consensuses[consensus_ID].nodes_IDs[consensus_nodes_count] = node_ID
                poagraph.consensuses[consensus_ID].nodes_count += 1

                #poagraph.consensuses[sequence_ID - max_source_ID - 1].add_node_ID(node_ID)

    for node_ID, line in enumerate(po_file_handler):
        base = line[0]
        in_nodes = _read_node_parameters(line, 'L')
        sequences_IDs = _read_node_parameters(line, 'S')
        aligned_to =_read_node_parameters(line, 'A')
        aligned_to = aligned_to[0] if aligned_to else None
        node = Node(ID=node_ID,
                    base=base,
                    in_nodes=np.array(in_nodes),
                    aligned_to=aligned_to
                    ) #todo bez consensus_count
        assign_this_node_to_its_sequences(sequences_IDs, node_ID)
        poagraph.add_node(node)

    return poagraph
#
# def _read_nodes_from_po_lines(poagraph, po_lines, max_source_ID):
#     def update_aligned_nodes_sets(aligned_node_sets, node_ID, node_aligned_to):
#         for s in aligned_node_sets:
#             if node_ID in s:
#                 s.add(node_aligned_to)
#         else:
#             aligned_node_sets.append(set([node_ID, node_aligned_to]))
#
#     def update_nodes_with_aligned_nodes(aligned_nodes_sets):
#         for s in aligned_nodes_sets:
#             for node_ID in s:
#                 poagraph.nodes[node_ID].aligned_to.update(s - set([node_ID]))
#
#     def assign_this_node_to_its_sequences(sequences_IDs, node_ID):
#         for sequence_ID in sequences_IDs:
#             if sequence_ID <= max_source_ID:
#                 poagraph.sources[sequence_ID].nodes_IDs.append(node_ID)
#             else:
#                 poagraph.consensuses[sequence_ID - max_source_ID - 1].add_node_ID(node_ID)
#
#     first_node_line_number = _get_first_node_line_number(po_lines)
#     aligned_nodes_sets = []
#     for i, line in enumerate(po_lines[first_node_line_number:]):
#         base = line[0]
#         in_nodes = set(_read_node_parameters(line, 'L'))
#         sequences_IDs = _read_node_parameters(line, 'S')
#         sources = set([sequence_ID for sequence_ID in sequences_IDs if sequence_ID <= max_source_ID])
#         consensuses_count = len(sequences_IDs) - len(sources)
#         poagraph.add_node(
#             Node(ID=i,
#                  base=base,
#                  in_nodes=in_nodes,
#                  sources = sources,
#                  consensuses_count=consensuses_count))
#
#         assign_this_node_to_its_sequences(sequences_IDs, i)
#
#         aligned_to = _read_node_parameters(line, 'A')
#         if aligned_to:
#             update_aligned_nodes_sets(aligned_nodes_sets, i, aligned_to[0])
#
#     update_nodes_with_aligned_nodes(aligned_nodes_sets)

# def read_single_consensus(file_path, consensusID=0):
#     print('\tRead consensus ' + str(consensusID) + ' from ' + file_path)
#     with open(file_path) as po:
#         po_lines = po.readlines()
#         p = POAGraph(name=_read_value(po_lines[1]),
#                      title=_read_value(po_lines[2]),
#                      version=_read_value(po_lines[0]),
#                      path='')
#         _read_sequence_info_from_po_lines(p, po_lines)
#         _read_nodes_from_po_lines(p, po_lines, len(p.sources)-1)
#     if not p.consensuses:
#         raise NoConsensusFound
#     else:
#         return p.consensuses[consensusID]
#
# def _spread_aligned_nodes(nodes):
#     column_alignment_cycle = set()
#     for n in nodes.values():
#         if not n.alignedTo:
#             continue
#         if min(n.alignedTo) in column_alignment_cycle:
#             n.alignedTo = set(column_alignment_cycle)
#             n.alignedTo.remove(n.ID)
#             continue
#         else:
#             aligned_node_id = max(n.alignedTo)
#             column_alignment_cycle = set([n.ID])
#             while True:
#                 column_alignment_cycle.add(aligned_node_id)
#                 next_aligned_node_ID = max(nodes[aligned_node_id].alignedTo)
#                 if next_aligned_node_ID == n.ID:
#                     break
#                 aligned_node_id = next_aligned_node_ID
#
#             n.alignedTo = set(column_alignment_cycle)
#             n.alignedTo.remove(n.ID)
#
#
# def _get_first_sequence_info_line_number():
#     return 5
#
#
# def _get_first_node_line_number(po_lines):
#     for i, l in enumerate(po_lines):
#         if l[1] == ':':
#             return (i)
#
#

