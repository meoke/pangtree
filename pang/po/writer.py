from graph.Pangraph import Pangraph
from pathlib import Path
from typing import List
from metadata import MultialignmentMetadata


def save(p: Pangraph, path: Path, genomes_info: MultialignmentMetadata) -> None:

    po_lines = [None] * get_poa_file_length(p)

    introduction = get_introduction(version=genomes_info.version,
                                    nodes_count=p.get_nodes_count(),
                                    source_count=p.get_paths_count(),
                                    title=genomes_info.title)
    file_cursor = len(introduction)
    po_lines[0:len(introduction)] = introduction

    sources = get_sources(p, genomes_info)
    po_lines[file_cursor: file_cursor + len(sources)] = sources
    file_cursor += len(sources)

    nodes = get_nodes(p)
    po_lines[file_cursor: file_cursor + len(nodes)] = nodes

    # for j, i in enumerate(po_lines):
    #     if i is None:
    #         po_lines[j] = ""

    # po_lines.extend(get_source_sequences(nodes, sources_weights))
    # n = get_nodes(nodes, sources)
    # po_lines.extend(n)

    lines = "\n".join(po_lines)
    with open(path, 'w') as out:
        out.write(lines)


def get_poa_file_length(pangraph: Pangraph) -> int:
    const_lines_count = 4
    source_lines_count = 4 * pangraph.get_paths_count()
    nodes_lines_count = pangraph.get_nodes_count()
    return const_lines_count + source_lines_count + nodes_lines_count


def get_introduction(version: str, nodes_count: int, source_count: int, title: str) -> List[str]:
    introduction_data = [f"VERSION={version}",
                         f"NAME={title}",
                         f"TITLE={title}",
                         f"LENGTH={str(nodes_count)}",
                         f"SOURCECOUNT={str(source_count)}"]
    return introduction_data


def get_sources(pangraph: Pangraph, genomes_metadata: MultialignmentMetadata) -> List[str]:
    source_count = pangraph.get_paths_count()
    sources_data = [None] * source_count * 2
    line_id = 0
    for source in pangraph.get_path_names():
        sources_data[line_id] = f"SOURCENAME={source}"
        sources_data[line_id + 1] = ("SOURCEINFO=" +
                                " ".join([f"{pangraph.get_path_nodes_count(source)}",
                                f"{pangraph.get_start_node_id(source)}",
                                f"{pangraph.get_source_weight(source)}",
                                f"{pangraph.get_source_consensus_id(source)}",
                                f"{source}"]))
        line_id += 2
    return sources_data


def get_nodes(pangraph: Pangraph) -> List[str]:
    for i in range(pangraph.get_nodes_count()):
