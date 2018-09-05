import re
import numpy as np
from POAGraph import POAGraph
from Sequence import Source, Consensus
from Node import Node
from Errors import NoConsensusFound


def parse_to_poagraph(file_path, output_dir):
    print('\tBuliding poagraph from ' + file_path) # todo logging
    with open(file_path) as po:
        po_lines = po.readlines()
    lines_iterator = iter(po_lines)
    poagraph = POAGraph(version=_extract_line_value(next(lines_iterator)),
                        name=_extract_line_value(next(lines_iterator)),
                        title=_extract_line_value(next(lines_iterator)),
                        path=output_dir)

    nodes_count = int(_extract_line_value(next(lines_iterator)))
    poagraph.sources, poagraph.consensuses = _read_sequence_info(lines_iterator)
    poagraph.nodes, poagraph.ns, poagraph.nc = _read_nodes_info(lines_iterator, nodes_count, len(poagraph.sources), len(poagraph.consensuses))
    return poagraph

def _extract_line_value(line):
    return line.split('=')[1].strip()

def _read_sequence_info(lines_iterator):
    source_count = int(_extract_line_value(next(lines_iterator)))

    source_ID = -1
    consensus_ID = -1
    sources = []
    consensuses = []
    for line in lines_iterator:
        sequence_name = _extract_line_value(line)
        detailed_info_line = next(lines_iterator)#po_file_handler.readline()
        detailed_info = _extract_line_value(detailed_info_line).split(' ')
        if 'CONSENS' in sequence_name:
            consensus_ID += 1
            consensus = Consensus(ID=consensus_ID,
                                  name=sequence_name,
                                  title=" ".join(detailed_info[4:]))
            consensuses.append(consensus)
        else:
            source_ID += 1
            source = Source(ID=source_ID,
                            name=sequence_name,
                            title=" ".join(detailed_info[4:]),
                            weight=int(detailed_info[2]))
            source.consensus_ID = int(detailed_info[3])
            sources.append(source)

        if source_ID + consensus_ID + 2 == source_count:
            break

    return sources, consensuses

def _read_nodes_info(lines_iterator, nodes_count, sources_count, consensuses_count):
    def assign_this_node_to_its_sequences(sequences_IDs, node_ID):
        sequeunces_IDs = np.array(sequences_IDs)
        srcs_IDs = sequeunces_IDs[sequeunces_IDs < sources_count]
        ns[srcs_IDs, node_ID] = True

        cons_ID = np.array([seq_ID - sources_count for seq_ID in sequeunces_IDs[sequeunces_IDs>=sources_count]])
        if cons_ID.size:
            nc[cons_ID, node_ID] = True


    nodes = [None] * nodes_count
    ns = np.zeros(shape=(sources_count, nodes_count), dtype=np.bool)
    nc = np.zeros(shape=(consensuses_count, nodes_count), dtype=np.bool)

    for node_ID, line in enumerate(lines_iterator):#po_file_handler):
        base = line[0].upper()
        in_nodes, sequences_IDs, aligned_to = _extract_node_parameters(line)
        # in_nodes = _extract_node_parameters(line, 'L')
        # sequences_IDs = _extract_node_parameters(line, 'S')
        # aligned_to =_extract_node_parameters(line, 'A')
        aligned_to = aligned_to[0] if aligned_to else None
        node = Node(ID=node_ID,
                    base=base,
                    in_nodes=np.array(in_nodes),
                    aligned_to=aligned_to
                    )
        assign_this_node_to_its_sequences(sequences_IDs, node_ID)
        nodes[node_ID] = node

    return nodes, ns, nc


def _extract_node_parameters(line):
    line = line[2:]
    line_iter = iter(line)

    node_parameters = {'L':[], 'S': [], 'A': [], '':[]}
    current_node_parameter = next(line_iter)
    number_start = 1
    number_end = 1
    for i, c in enumerate(line_iter):
        if c == 'L' or c == 'S' or c == 'A':
            node_parameters[current_node_parameter].append(int(line[number_start: number_end]))
            current_node_parameter = c
            number_start = i + 2
            number_end = i + 2
        else:
            number_end += 1
    node_parameters[current_node_parameter].append(int(line[number_start: number_end]))
    return node_parameters['L'], node_parameters['S'], node_parameters['A']
    # pattern = '{0}\d+'.format(code_letter)
    # values_with_prefix_letters = re.findall(pattern, node)
    # return [int(letter_value[1:]) for letter_value in values_with_prefix_letters]


def read_consensus(po_file_path, consensusID=0):
    print('\tRead consensus ' + str(consensusID) + ' from ' + po_file_path)
    poagraph = parse_to_poagraph(po_file_path, output_dir="")
    if not poagraph.consensuses:
        raise NoConsensusFound
    else:
        return poagraph.consensuses[consensusID], poagraph.nc[consensusID][:]
