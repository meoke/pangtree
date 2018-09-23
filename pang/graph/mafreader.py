from pathlib import Path
from Bio import AlignIO
import numpy as np
from .Graph import Graph
from .PathManager import PathManager
from maf.Mafgraph import Mafgraph
from metadata.MultialignmentMetadata import MultialignmentMetadata
from .Node import Node
from graph import nucleotides

# marfreader - jako że maf zawiera też info o ścieżkach, no to może zwracać też pathmanagera!

def read(multialignment_file: Path, multialignment_metadata: MultialignmentMetadata):
    maf = AlignIO.parse(multialignment_file, "maf")
    return parse(maf, multialignment_metadata)


def parse(maf: str, multialignment_metadata):
    mafgraph = get_sorted_mafgraph(maf) #żeby znać kolejność
    nodes_count = get_max_count_nodes(mafgraph)
    sequences_count = len(multialignment_metadata.genomes_metadata)

    graph = Graph(nodes_count=nodes_count)
    paths_names=[m.mafname for m in multialignment_metadata.genomes_metadata.values()]
    paths_manager = PathManager(start_node_id=0, max_nodes_count=nodes_count, paths_names=paths_names)
    last_node_id = -1
    for b in mafgraph:
        nodes, minipathmanager, last_node_id = parse_maf_block(b, last_node_id)
        graph.update_nodes(nodes)
        paths_manager.update(minipathmanager, last_node_id+1-len(nodes))
    pass


def get_sorted_mafgraph(maf):
    return Mafgraph(maf=maf, remove_cycles=True)


def get_max_count_nodes(mafgraph):
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
    block_sequences = [seq.id for seq in block.alignment]
    minipathmanager = PathManager(start_node_id=previous_node_id+1,
                                  max_nodes_count=block_width * 4,
                                  paths_names=block_sequences)

    nodes = []
    current_node_id = previous_node_id
    for col in range(block_width):
        d = {seq.id: seq[col] for seq in block.alignment}
        nodes_codes = [*(set(d.values())).difference(set(['-']))]

        for j, n in enumerate(nodes_codes):
            current_node_id += 1
            node = Node(id=current_node_id,
                        base=nucleotides.code(n),
                        in_nodes=[],
                        aligned_to=get_next_aligned_node_id(j, nodes_codes, previous_node_id))
            for sequence, nucleotide in d.items():
                if nucleotide == n:
                    minipathmanager.mark(path_name=sequence, node_id=current_node_id)
            nodes.append(node)

    minipathmanager.trim(current_node_id)
    return nodes, minipathmanager, current_node_id