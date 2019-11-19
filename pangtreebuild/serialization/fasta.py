from pangtreebuild.affinity_tree import tree
from pangtreebuild.pangenome import graph


def poagraph_to_fasta(poagraph: graph.Poagraph) -> str:
    """Converts poagraph to FASTA format.

    Args:
        poagraph: Poagraph to be converted.

    Returns:
        Fasta formatted string ready to be saved to file.
    """

    fasta_lines = []
    for seq_id, sequence in poagraph.sequences.items():
        if poagraph.get_sequence_nodes_count(seq_id) == 0:
            continue
        sequence = "".join([poagraph.nodes[node_id].get_base()
                            for path in sequence.paths
                            for node_id in path])
        fasta_lines.append(f">{seq_id}")
        fasta_lines.append(sequence)

    return "\n".join(fasta_lines)


def affinity_tree_to_fasta(poagraph: graph.Poagraph,
                           affinity_tree: tree.AffinityTree) -> str:
    """Converts Affinity tree to FASTA formatted string.

    Args:
        poagraph: Poagraph linked to the affinity_tree.
        affinity_tree: Consensus paths from this affinty tree
            will be converted to FASTA.

    Returns:
        Fasta formatted string with consensus paths from given affinty tree.
    """

    fasta_lines = []
    for affinity_node in affinity_tree.nodes:
        if affinity_node.consensus is None or \
           len(affinity_node.consensus) == 0:
            continue
        sequence = "".join([poagraph.nodes[node_id].get_base()
                            for node_id in affinity_node.consensus])
        missing_mincomp_symbol = "?"
        leaf = "" \
               if len(affinity_node.children) > 1 \
               else ",".join([str(seq_id)
                              for seq_id in affinity_node.sequences])
        fasta_lines.append(f">AffinityNode{affinity_node.id_}|"
                           f"leaf_for={leaf}|"
                           f"""mincomp={affinity_node.mincomp
                                      if affinity_node.mincomp is not None
                                      else missing_mincomp_symbol}|"""
                           f"sequences_count={len(affinity_node.sequences)}|"
                           f"children={str(affinity_node.children)}")
        fasta_lines.append(sequence)

    return "\n".join(fasta_lines)
