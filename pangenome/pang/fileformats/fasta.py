from consensus.ConsensusesTree import ConsensusesTree
from pangraph.Pangraph import Pangraph


def pangraph_to_fasta(pangraph: Pangraph) -> str:
    fasta_lines = []
    for seq_id, paths in pangraph.paths.items():
        if pangraph.path_is_empty(seq_id):
            continue
        sequence = "".join([pangraph.nodes[node_id].get_base() for path in paths for node_id in path])
        fasta_lines.append(f">{seq_id}")
        fasta_lines.append(sequence)

    return "\n".join(fasta_lines)


def consensuses_tree_to_fasta(pangraph: Pangraph, consensuses_tree: ConsensusesTree) -> str:
    fasta_lines = []
    for consensus_node in consensuses_tree.nodes:
        if len(consensus_node.consensus_path) == 0:
            continue
        sequence = "".join([pangraph.nodes[node_id].get_base() for node_id in consensus_node.consensus_path])
        missing_mincomp_symbol = "?"
        fasta_lines.append(f">consensus{consensus_node.consensus_id}|"
                           f"mincomp={consensus_node.mincomp if consensus_node.mincomp is not None else missing_mincomp_symbol}|"
                           f"sequences_count={len(consensus_node.sequences_ids)}|"
                           f"children={str(consensus_node.children_nodes_ids)}")
        fasta_lines.append(sequence)

    return "\n".join(fasta_lines)
