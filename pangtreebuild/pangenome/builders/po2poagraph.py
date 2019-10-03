from collections import namedtuple
from typing import List, Optional, Tuple, Dict

from pangtreebuild.pangenome import poagraph
from pangtreebuild.pangenome.parameters import multialignment

# global_logger = loggingtools.get_global_logger()
# detailed_logger = loggingtools.get_logger("details")

POSequenceInfo = namedtuple('POSequenceInfo', ['name',
                                               'nodes_count',
                                               'start_node',
                                               'weight',
                                               'consensus_id',
                                               'additional_info'])


def get_poagraph(po: multialignment.Po, metadataCsv: multialignment.MetadataCSV):
    po_lines = po.filecontent.readlines()
    sequences_info = _get_sequences_info_from_po(po_lines)
    initial_sequences = _init_sequences(sequences_info, metadataCsv)
    nodes, sequences = _get_poagraph_paths_and_nodes(po_lines, sequences_info, initial_sequences)
    return nodes, sequences


def _extract_line_value(line: str) -> str:
    splitted_line = line.split('=')
    if len(splitted_line) < 2:
        raise Exception("Incorrect line format. Nothing after \"=\" found.")
    return splitted_line[1].strip()


def _init_sequences(sequences_info: Dict[int, POSequenceInfo],
                    metadata: Optional[multialignment.MetadataCSV]) -> Dict[multialignment.SequenceID, poagraph.Sequence]:
    metadata_sequences_ids = metadata.get_all_sequences_ids() if metadata else []
    po_sequences_ids = [seq_info.name for seq_info in sequences_info.values()]
    initial_sequences = {seq_id: poagraph.Sequence(seqid=seq_id,
                                          paths=[],
                                          seqmetadata=metadata.get_sequence_metadata(seq_id)
                                          if metadata else {})
                         for seq_id in set(po_sequences_ids + metadata_sequences_ids)}

    return initial_sequences


def _get_sequences_info_from_po(po_lines: List[str]) -> Dict[int, POSequenceInfo]:
    paths_count = int(_extract_line_value(po_lines[4]))
    sequences_info = {}
    po_seq_id = 0
    for i in range(5, 5 + paths_count * 2, 2):
        path_name = _extract_line_value(po_lines[i])
        detailed_info_line = po_lines[i+1]
        detailed_info = _extract_line_value(detailed_info_line).split(' ')
        if len(detailed_info) < 5:
            raise Exception(f"Expeceted SOURCEINFO=[5 parameters]. Got {detailed_info_line} instead.")
        sequences_info[po_seq_id] = POSequenceInfo(name=multialignment.SequenceID(path_name),
                                                   nodes_count=detailed_info[0],
                                                   start_node=detailed_info[1],
                                                   weight=detailed_info[2],
                                                   consensus_id=detailed_info[3],
                                                   additional_info="".join(detailed_info[4:]))
        po_seq_id += 1
    return sequences_info


def _get_poagraph_paths_and_nodes(po_lines: List[str],
                                  sequences_info: Dict[int, POSequenceInfo],
                                  sequences: Dict[multialignment.SequenceID, poagraph.Sequence]) -> \
        Tuple[List[poagraph.Node], Dict[multialignment.SequenceID, poagraph.Sequence]]:
    nodes_count = int(_extract_line_value(po_lines[3]))
    paths_count = int(_extract_line_value(po_lines[4]))
    #todo uzupelniac kolumn id_
    nodes: List[poagraph.Node] = [None] * nodes_count
    node_id = 0
    for i in range(5 + paths_count * 2, 5 + paths_count * 2 + nodes_count):
        node_line = po_lines[i]
        base = poagraph.Base(node_line[0].upper())
        in_nodes, po_sequences_ids, aligned_to = _extract_node_parameters(node_line)
        sequences_ids = [sequences_info[po_sequences_id].name for po_sequences_id in po_sequences_ids]
        nodes[node_id] = poagraph.Node(poagraph.NodeID(node_id), base, poagraph.NodeID(aligned_to))
        for seq_id in sequences_ids:
            if len(sequences[seq_id].paths) == 1:
                sequences[seq_id].paths[0].append(poagraph.NodeID(node_id))
            else:
                sequences[seq_id].paths.append(poagraph.SeqPath([poagraph.NodeID(node_id)]))
        node_id += 1
    return nodes, sequences


def _extract_node_parameters(line: str) -> Tuple[List[int], List[int], int]:
    line = line[2:]
    line_iter = iter(line)

    node_parameters = {'L': [], 'S': [], 'A': [], '': []}
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
