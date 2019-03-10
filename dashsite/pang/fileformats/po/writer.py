from pangraph.Pangraph import Pangraph
from pathlib import Path
from typing import List
from metadata.MultialignmentMetadata import MultialignmentMetadata
from pangraph.Node import Node
from pangraph import nucleotides as n


def save(p: Pangraph, path: Path, genomes_info: MultialignmentMetadata) -> None:
    raise NotImplementedError()

    po_lines = [None] * get_poa_file_length(p)

    introduction = get_introduction(version=genomes_info.version,
                                    nodes_count=p.get_nodes_count(),
                                    source_count=p.get_paths_count(),
                                    title=genomes_info.title)
    file_cursor = len(introduction)
    po_lines[0:len(introduction)] = introduction

    sources = get_sources(p)
    po_lines[file_cursor: file_cursor + len(sources)] = sources
    file_cursor += len(sources)

    nodes = get_nodes(p)
    po_lines[file_cursor: file_cursor + len(nodes)] = nodes

    lines = "\n".join(po_lines)
    with open(path, 'w') as out:
        out.write(lines)


def get_poa_file_length(pangraph: Pangraph) -> int:
    const_lines_count = 4
    source_lines_count = 2 * pangraph.get_paths_count()
    nodes_lines_count = pangraph.get_nodes_count()
    return const_lines_count + source_lines_count + nodes_lines_count


def get_introduction(version: str, nodes_count: int, source_count: int, title: str) -> List[str]:
    introduction_data = [f"VERSION={version}",
                         f"NAME={title}",
                         f"TITLE={title}",
                         f"LENGTH={str(nodes_count)}",
                         f"SOURCECOUNT={str(source_count)}"]
    return introduction_data


def get_sources(pangraph: Pangraph) -> List[str]:
    source_count = pangraph.get_paths_count()
    sources_data = [None] * source_count * 2
    line_id = 0
    sources_weights = pangraph.get_sources_weights_dict()
    for source_name, weight in sources_weights.items():
        sources_data[line_id] = f"SOURCENAME={source_name}"
        sources_data[line_id + 1] = ("SOURCEINFO=" +
                                     " ".join([f"{pangraph.get_path_nodes_count(source_name)}",
                                               f"{pangraph.get_start_node_id(source_name)}",
                                               f"{weight}",
                                               f"{pangraph.get_source_consensus_id(source_name)}",
                                               f"{source_name}"]))
        line_id += 2
    return sources_data


def get_node_code(node):
    return n.decode(node.base).lower()


def get_nodes(pangraph: Pangraph) -> List[str]:
    nodes_count = pangraph.get_nodes_count()
    nodes_data = [None] * nodes_count
    # todo write as list comprehension
    for i, node in enumerate(pangraph.get_nodes()):
        sources_ids = pangraph.get_sources_ids(node.id)
        nodes_data[i] = "".join([get_node_code(node),
                                 ":",
                                 get_in_nodes_info(node),
                                 get_sources_info(sources_ids),
                                 get_aligned_to_info(node)])
    return nodes_data


def get_in_nodes_info(node: Node) -> str:
    return "".join([f'L{i}' for i in node.in_nodes])


def get_sources_info(sources_ids):
    return "".join([f'S{i}' for i in sources_ids])


def get_aligned_to_info(node):
    return f"A{node.aligned_to}" if node.aligned_to is not None else ""
