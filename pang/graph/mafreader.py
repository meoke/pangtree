from pathlib import Path
from Bio import AlignIO
from .Graph import Graph
from .PathManager import PathManager
from maf.Mafgraph import Mafgraph
from metadata.MultialignmentMetadata import MultialignmentMetadata
from .Node import Node
from graph import nucleotides


def read(multialignment_file: Path, multialignment_metadata: MultialignmentMetadata) -> (Graph, PathManager):
    maf = AlignIO.parse(multialignment_file, "maf")
    return parse_mafcontent_to_graph(maf, multialignment_metadata)


def parse_mafcontent_to_graph(maf: str, multialignment_metadata: MultialignmentMetadata) -> (Graph, PathManager):
    mafgraph = get_sorted_mafgraph(maf)
    max_nodes_count = get_max_count_nodes(mafgraph)
    path_manager = PathManager(start_node_id=0,
                               max_nodes_count=max_nodes_count,
                               paths_names=[seq.mafname for seq in multialignment_metadata.genomes_metadata.values()])

    graph = Graph(nodes_count=max_nodes_count)
    last_node_id = -1
    for mafblock in mafgraph:
        nodes, block_path_manager, last_node_id = parse_maf_block(mafblock, last_node_id)
        graph.update_nodes(nodes)
        path_manager.update(block_path_manager, last_node_id+1-len(nodes))
    graph.trim(nodes_count=last_node_id+1)
    path_manager.trim(nodes_count=last_node_id+1)
    in_nodes_setup(graph, path_manager)
    return graph, path_manager

def in_nodes_setup(graph, path_manager):
    for node in graph.nodes:
        node.in_nodes = path_manager.get_in_nodes(node.id)

def get_sorted_mafgraph(maf) -> Mafgraph:
    #todo typ maf (generator?)
    return Mafgraph(maf=maf, remove_cycles=True)


def get_max_count_nodes(mafgraph: Mafgraph) -> int:
    #todo czy może być takie tymczasowe rozwiązanie:
    #szerokość każdego bloku maf *4 (bo to maksymalna liczba węzłów) i potem przyciąć tablicę

    max_nodes_count = 0
    for block in mafgraph.blocks:
        block_width = len(block.alignment[0].seq)
        max_nodes_count += block_width * 4
    return max_nodes_count


def get_next_aligned_node_id(j, nodes_codes, previous_node_id):
    if len(nodes_codes) > 1:
        return previous_node_id + 1 + ((j + 1) % len(nodes_codes))
    return None


def parse_maf_block(block, previous_node_id):
    block_width = len(block.alignment[0].seq)
    block_sequences_names = [seq.id for seq in block.alignment]
    block_pm = PathManager(start_node_id=previous_node_id+1,
                           max_nodes_count=block_width * 4,
                           paths_names=block_sequences_names)

    current_node_id = previous_node_id
    nodes = []
    for col in range(block_width):
        sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block.alignment}
        nodes_codes = [*(set(sequence_name_to_nucleotide.values())).difference(set(['-']))]

        for i, nucl in enumerate(nodes_codes):
            current_node_id += 1
            node = Node(id=current_node_id,
                        base=nucleotides.code(nucl),
                        in_nodes=[],
                        aligned_to=get_next_aligned_node_id(i, nodes_codes, previous_node_id))
            for sequence, nucleotide in sequence_name_to_nucleotide.items():
                if nucleotide == nucl:
                    block_pm.mark(path_name=sequence, node_id=current_node_id)
            nodes.append(node)



    block_pm.trim(nodes_count=len(nodes))

    return nodes, block_pm, current_node_id
