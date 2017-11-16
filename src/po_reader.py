import re
from POAGraph import POAGraph
from Sequence import Consensus, Source
from Node import Node
from NoConsensusFound import *

def parse_to_poagraph(file_path, output_dir):
    print('\tBuliding poagraph from ' + file_path)
    with open(file_path) as po:
        po_lines = po.readlines()
        p = POAGraph(name = _read_value(po_lines[1]),
                     title = _read_value(po_lines[2]),
                     version = _read_value(po_lines[0]),
                     path = output_dir)
        _read_sequence_info_from_po_lines(p, po_lines)
        _read_nodes_from_po_lines(p, po_lines, len(p.sources)-1)
    return p

def read_single_consensus(file_path, consensusID=0):
    print('\tRead consensus ' + str(consensusID) + ' from ' + file_path)
    with open(file_path) as po:
        po_lines = po.readlines()
        p = POAGraph(name = _read_value(po_lines[1]),
                     title = _read_value(po_lines[2]),
                     version = _read_value(po_lines[0]),
                     path = '')
        _read_sequence_info_from_po_lines(p, po_lines)
        _read_nodes_from_po_lines(p, po_lines, len(p.sources)-1)
    if not p.consensuses:
        raise NoConsensusFound
    else:
        return p.consensuses[consensusID]

def _read_value(line):
    return line.split('=')[1].strip()


def _read_sequence_info_from_po_lines(poagraph, po_lines):
    sequence_info_start = _get_first_sequence_info_line_number()
    node_info_start = _get_first_node_line_number(po_lines)
    source_ID = -1
    consensus_ID = -1
    lines_iterator = iter(po_lines[sequence_info_start:node_info_start])
    for line in lines_iterator:
        sequence_name = _read_value(line)

        if 'CONSENS' in sequence_name:
            consensus_ID += 1
            consensus = Consensus(currentID=consensus_ID, name=sequence_name, title=sequence_name)
            line = next(lines_iterator)
            consensus_info = _read_value(line).split(' ')
            consensus.title = " ".join(consensus_info[4:])
            poagraph.add_consensus(consensus)
        else:
            source_ID += 1
            source = Source(currentID=source_ID, name=sequence_name, title=sequence_name)
            line = next(lines_iterator)
            source_info = _read_value(line).split(' ')
            source.weight = int(source_info[2])
            source.consensusID = int(source_info[3])
            source.title = " ".join(source_info[4:])
            poagraph.add_source(source)


def _read_nodes_from_po_lines(poagraph, po_lines, max_source_ID):
    def update_aligned_nodes_sets(aligned_node_sets, node_ID, node_aligned_to):
        for s in aligned_node_sets:
            if node_ID in s:
                s.add(node_aligned_to)
        else:
            aligned_node_sets.append(set([node_ID, node_aligned_to]))

    def update_nodes_with_aligned_nodes(aligned_nodes_sets):
        for s in aligned_nodes_sets:
            for node_ID in s:
                poagraph.nodes[node_ID].aligned_to.update(s - set([node_ID]))

    def assign_this_node_to_its_sequences(sequences_IDs, node_ID):
        for sequence_ID in sequences_IDs:
            if sequence_ID <= max_source_ID:
                poagraph.sources[sequence_ID].nodes_IDs.append(node_ID)
            else:
                poagraph.consensuses[sequence_ID - max_source_ID - 1].add_node_ID(node_ID)

    first_node_line_number = _get_first_node_line_number(po_lines)
    aligned_nodes_sets = []
    for i, line in enumerate(po_lines[first_node_line_number:]):
        base = line[0]
        in_nodes = set(_read_node_parameters(line, 'L'))
        sequences_IDs = _read_node_parameters(line, 'S')
        sources = set([sequence_ID for sequence_ID in sequences_IDs if sequence_ID <= max_source_ID])
        consensuses_count = len(sequences_IDs) - len(sources)
        poagraph.add_node(
            Node(currentID=i,
            base=base,
            in_nodes=in_nodes,
            sources = sources,
            consensuses_count=consensuses_count))

        assign_this_node_to_its_sequences(sequences_IDs, i)

        aligned_to = _read_node_parameters(line, 'A')
        if aligned_to:
            update_aligned_nodes_sets(aligned_nodes_sets, i, aligned_to[0])

    update_nodes_with_aligned_nodes(aligned_nodes_sets)

def _spread_aligned_nodes(nodes):
    column_alignment_cycle = set()
    for n in nodes.values():
        if not n.alignedTo:
            continue
        if min(n.alignedTo) in column_alignment_cycle:
            n.alignedTo = set(column_alignment_cycle)
            n.alignedTo.remove(n.ID)
            continue
        else:
            aligned_node_id = max(n.alignedTo)
            column_alignment_cycle = set([n.ID])
            while True:
                column_alignment_cycle.add(aligned_node_id)
                next_aligned_node_ID = max(nodes[aligned_node_id].alignedTo)
                if next_aligned_node_ID == n.ID:
                    break
                aligned_node_id = next_aligned_node_ID

            n.alignedTo = set(column_alignment_cycle)
            n.alignedTo.remove(n.ID)


def _get_first_sequence_info_line_number():
    return 5


def _get_first_node_line_number(po_lines):
    for i, l in enumerate(po_lines):
        if l[1] == ':':
            return (i)


def _read_node_parameters(node, code_letter):
    pattern = '{0}\d+'.format(code_letter)
    values_with_prefix_letters = re.findall(pattern, node)
    return [int(letter_value[1:]) for letter_value in values_with_prefix_letters]
