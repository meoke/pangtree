from typing import List, Tuple, Callable

from pangtreebuild.pangenome import graph


class NodePO:
    """Poagraph node representation in po file."""

    def __init__(self,
                 base: bytes,
                 aligned_to: graph.NodeID,
                 in_nodes: List[graph.NodeID],
                 sequences_ids: List[int]):
        self.base = base
        self.aligned_to = aligned_to
        self.in_nodes = in_nodes
        self.sequences_ids = sequences_ids


class SequencePO:
    """Sequences representation in po file."""

    def __init__(self,
                 name: str,
                 nodes_count: int,
                 weight: int,
                 consensus_id: int,
                 start_node_id: graph.NodeID):
        self.name = name
        self.nodes_count = nodes_count
        self.weight = weight
        self.consensus_id = consensus_id
        self.start_node_id = start_node_id


def poagraph_to_PangenomePO(poagraph: graph.Poagraph) -> str:
    """Converts poagraph to .po file."""

    po_nodes, po_sequences = _convert_to_po_input_data(poagraph)
    return poagraph_elements_to_PangenomePO(po_nodes,
                                            po_sequences,
                                            poagraph.datatype)


def poagraph_elements_to_PangenomePO(po_nodes: List[NodePO],
                                     po_sequences: List[SequencePO],
                                     datatype: graph.DataType) -> str:
    """Converts po specific elements to po file."""

    introduction_lines = "\n".join(_get_introduction_lines(len(po_nodes),
                                                           len(po_sequences)))
    sequences_lines = "\n".join(_get_sequences_lines(po_sequences))
    nodes_lines = "\n".join(_get_nodes_lines(po_nodes, datatype))
    return "\n".join([introduction_lines, sequences_lines, nodes_lines])


def _convert_to_po_input_data(p: graph.Poagraph) -> \
        Tuple[List[NodePO], List[SequencePO]]:
    po_nodes = []
    po_sequences = []
    sequences_weights = p.get_sequences_weights(p.get_sequences_ids())

    for node in p.nodes:
        po_nodes.append(NodePO(base=node._base.value,
                               aligned_to=node.aligned_to,
                               in_nodes=set(),
                               sequences_ids=[]))

    seq_int_id = -1
    for seq_id, sequence in p.sequences.items():
        nodes_count = sum([len(path) for path in sequence.paths])
        if nodes_count == 0:
            continue
        seq_int_id += 1
        po_sequences.append(SequencePO(name=str(seq_id),
                                       nodes_count=nodes_count,
                                       weight=sequences_weights[seq_id],
                                       consensus_id=-1,
                                       start_node_id=sequence.paths[0][0]))
        for path in sequence.paths:
            previous_node_id = None
            for node_id in path:
                po_nodes[node_id].sequences_ids.append(seq_int_id)
                if previous_node_id is not None:
                    po_nodes[node_id].in_nodes.add(previous_node_id)
                previous_node_id = node_id

    for node in po_nodes:
        node.in_nodes = list(node.in_nodes)

    return po_nodes, po_sequences


def _get_introduction_lines(nodes_count: int, sequences_count: int) -> \
        List[str]:
    introduction_lines = ["VERSION=pangenome",
                          "NAME=pangenome",
                          "TITLE=pangenome",
                          f"LENGTH={nodes_count}",
                          f"SOURCECOUNT={sequences_count}"]
    return introduction_lines


def _get_sequences_lines(po_sequences: List[SequencePO]) -> List[str]:
    if not po_sequences:
        raise Exception("No _sequences info to write in PO file.")
    i = -2
    po_sequences_lines = [None] * len(po_sequences) * 2
    for sequence in po_sequences:
        i += 2
        po_sequences_lines[i] = f"SOURCENAME={sequence.name}"
        po_sequences_lines[i + 1] = ("SOURCEINFO=" +
                                     " ".join([f"{sequence.nodes_count}",
                                               f"{sequence.start_node_id}",
                                               f"{sequence.weight}",
                                               f"{sequence.consensus_id}",
                                               f"{sequence.name}"]))
    return po_sequences_lines


def _get_nodes_lines(po_nodes: List[NodePO], datatype: graph.DataType) -> \
        List[str]:
    if not po_nodes:
        raise Exception("No nodes info to write in PO file.")
    i = 0

    _get_node_code = _get_node_code_conversion_function(datatype)
    po_nodes_lines = [None] * len(po_nodes)
    for node in po_nodes:
        po_nodes_lines[i] = "".join([_get_node_code(node.base),
                                     ":",
                                     _get_in_nodes_info(node.in_nodes),
                                     _get_sources_info(node.sequences_ids),
                                     _get_aligned_to_info(node.aligned_to)])
        i += 1
    return po_nodes_lines


def _get_in_nodes_info(in_nodes: List[graph.NodeID]) -> str:
    return "".join([f'L{i}' for i in in_nodes])


def _get_sources_info(sources_ids: List[int]) -> str:
    return "".join([f'S{i}' for i in sources_ids])


def _get_aligned_to_info(aligned_to: graph.NodeID) -> str:
    return f"A{aligned_to}" if aligned_to is not None else ""


def _get_node_code_conversion_function(datatype: graph.DataType) -> \
        Callable[[bytes], str]:
    if datatype.value == graph.DataType.Proteins.value:
        return _get_protein_node_code
    elif datatype.value == graph.DataType.Nucleotides.value:
        return _get_nucleotides_node_code
    else:
        raise Exception("Unknown data type. Cannot create PO file.")


def _get_protein_node_code(base: bytes) -> str:
    return base.decode("ASCII").upper()


def _get_nucleotides_node_code(base: bytes) -> str:
    return base.decode("ASCII").lower()
