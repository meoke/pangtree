from pathlib import Path
from Bio import AlignIO
from .Graph import Graph
from .PathManager import PathManager
from maf.Mafgraph import Mafgraph
from metadata.MultialignmentMetadata import MultialignmentMetadata
from .Node import Node
from graph import nucleotides
from .Pangraph import Pangraph


def read(multialignment_file: Path, multialignment_metadata: MultialignmentMetadata) -> Pangraph:
    maf = AlignIO.parse(multialignment_file, "maf")
    return parse_mafcontent_to_graph(maf, multialignment_metadata)


def parse_mafcontent_to_graph(maf: str, metadata: MultialignmentMetadata) -> Pangraph:
    mafgraph = get_sorted_mafgraph(maf)
    max_nodes_count = get_max_count_nodes(mafgraph)
    pangraph = Pangraph(max_nodes_count=max_nodes_count,
                        start_node_id=0,
                        paths_names=[seq.mafname for seq in metadata.genomes_metadata.values()])

    last_node_id = -1
    for mafblock in mafgraph:
        block_pangraph, last_node_id = parse_maf_block(mafblock, last_node_id)
        pangraph.update(block_pangraph, last_node_id+1-block_pangraph.get_nodes_count())
    pangraph.trim_nodes(nodes_count=last_node_id+1)
    pangraph.fill_in_nodes()
    pangraph.remove_empty_paths()
    return pangraph


def fill_in_nodes(graph: Graph, path_manager: PathManager):
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


def get_next_aligned_node_id(current_column_i, column_nodes_ids):
    if len(column_nodes_ids) > 1:
        return column_nodes_ids[(current_column_i+1) % len(column_nodes_ids)]
    return None


def parse_maf_block(block, previous_node_id):
    # todo uwzględniać + i - z sekwencji
    # todo uwzględniać czy blok jest odwrócony czy nie

    block_width = len(block.alignment[0].seq)
    block_sequences_names = [seq.id for seq in block.alignment]

    pg = Pangraph(max_nodes_count=block_width*4,
                  start_node_id=previous_node_id+1,
                  paths_names=block_sequences_names)
    current_node_id = previous_node_id
    for col in range(block_width):
        sequence_name_to_nucleotide = {seq.id: seq[col] for seq in block.alignment}
        nodes_codes = sorted([*(set(sequence_name_to_nucleotide.values())).difference(set(['-']))])
        column_nodes_ids = [current_node_id + i + 1 for i, _ in enumerate(nodes_codes)]

        for i, nucl in enumerate(nodes_codes):
            current_node_id += 1
            node = Node(id=current_node_id,
                        base=nucleotides.code(nucl),
                        in_nodes=[],
                        aligned_to=get_next_aligned_node_id(i, column_nodes_ids))
            for sequence, nucleotide in sequence_name_to_nucleotide.items():
                if nucleotide == nucl:
                    pg.add_path_to_node(path_name=sequence, node_id=current_node_id)
            pg.add_node(node, current_node_id-previous_node_id-1)
    pg.trim_nodes(nodes_count=current_node_id-previous_node_id)
    return pg, current_node_id
