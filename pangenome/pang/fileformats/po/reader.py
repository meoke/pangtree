from pathlib import Path

from pangraph.Pangraph import Pangraph
from pangraph.Pangraph import Node


def read(path: Path) -> Pangraph:
    raise NotImplementedError()

    with open(path) as po_input:
        po_lines = po_input.readlines()

    p = Pangraph(pangenome_parameters.datatype)
    nodes, sequences_to_nodes_ids, consensuses_to_nodes_ids = extract_pangraph(po_lines)
    p.update_nodes(new_nodes=nodes)
    p.set_paths(max_nodes_count=len(nodes), paths_to_node_ids=sequences_to_nodes_ids)
    p.set_consensuses(max_nodes_count=len(nodes), paths_to_node_ids=consensuses_to_nodes_ids)
    return p


def extract_pangraph(po_lines):
    po_lines_iterator = iter(po_lines)

    for i in range(3):
        next(po_lines_iterator)

    nodes_count = int(_extract_line_value(next(po_lines_iterator)))
    paths_count = int(_extract_line_value(next(po_lines_iterator)))

    seqnames_to_nodes_ids = {}
    consensusnames_to_nodes_ids = {}

    path_ids_to_pathnames = {}

    for i in range(paths_count):
        path_name = _extract_line_value(next(po_lines_iterator))
        path_ids_to_pathnames[i] = path_name
        detailed_path_info = next(po_lines_iterator)
        detailed_info = _extract_line_value(detailed_path_info).split(' ')

        path_nodes_count = int(detailed_info[0])
        consensus_number = int(detailed_info[3]) #todo use this info

        if 'CONSENS' in path_name:
            consensusnames_to_nodes_ids[path_name] = []# todo would be faster? [None] * path_nodes_count
        else:
            seqnames_to_nodes_ids[path_name] = []# todo would be faster? [None] * path_nodes_count

    nodes = [None] * nodes_count
    for i in range(nodes_count):
        node_line = next(po_lines_iterator)
        base = node_line[0].upper()
        in_nodes, sequences_IDs, aligned_to = _extract_node_parameters(node_line)
        nodes[i] = Node(node_id=i,
                        base=n.code(base),
                        in_nodes=in_nodes,
                        aligned_to=aligned_to,
                        column_id=-1,  #TODO
                        block_id=0)
        for sequence_id in sequences_IDs:
            pathname = path_ids_to_pathnames[sequence_id]
            if sequence_id < len(seqnames_to_nodes_ids.keys()):
                seqnames_to_nodes_ids[pathname].append(i)
            else:
                consensusnames_to_nodes_ids[pathname].append(i)

    return nodes, seqnames_to_nodes_ids, consensusnames_to_nodes_ids


def _extract_line_value(line):
    return line.split('=')[1].strip()


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
    aligned_to = node_parameters['A'][0] if node_parameters['A'] else None
    return node_parameters['L'], node_parameters['S'], aligned_to
