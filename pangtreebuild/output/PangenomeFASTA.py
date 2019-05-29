from poapangenome.consensus.ConsensusTree import ConsensusTree
from poapangenome.datamodel.Poagraph import Poagraph


def poagraph_to_fasta(poagraph: Poagraph) -> str:
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


def consensuses_tree_to_fasta(poagraph: Poagraph, consensus_tree: ConsensusTree) -> str:
    fasta_lines = []
    for consensus_node in consensus_tree.nodes:
        if consensus_node.consensus_path is None or len(consensus_node.consensus_path) == 0:
            continue
        sequence = "".join([poagraph.nodes[node_id].get_base()
                            for node_id in consensus_node.consensus_path])
        missing_mincomp_symbol = "?"
        fasta_lines.append(f">consensus{consensus_node.consensus_id}|"
                           f"mincomp={consensus_node.mincomp if consensus_node.mincomp is not None else missing_mincomp_symbol}|"
                           f"sequences_count={len(consensus_node.sequences_ids)}|"
                           f"children={str(consensus_node.children_nodes_ids)}")
        fasta_lines.append(sequence)

    return "\n".join(fasta_lines)
