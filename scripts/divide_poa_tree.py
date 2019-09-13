from pathlib import Path
from typing import Tuple

from build.lib.pangtreebuild.datamodel.Node import Node
from pangtreebuild.consensus.ConsensusTree import ConsensusTree, ConsensusNode, ConsensusNodeID, CompatibilityToPath
from pangtreebuild.datamodel.Node import NodeID, Base
from pangtreebuild.datamodel.Poagraph import Poagraph
from pangtreebuild.datamodel.Sequence import SequenceID, Sequence
from pangtreebuild.output.PangenomeJSON import str_to_PangenomeJSON
from pangtreebuild.tools import pathtools


def convert_jsonpangenome(path: Path) -> Tuple[Poagraph, ConsensusTree]:
    jsonpangenome = str_to_PangenomeJSON(pathtools.get_file_content(path))
    poagraph = Poagraph(nodes=[Node(node_id=NodeID(n.id), base=Base(n.base)) for n in jsonpangenome.nodes],
                    sequences={SequenceID(seq.sequence_str_id): Sequence(seqid=SequenceID(seq.sequence_str_id), paths=[], seqmetadata=seq.metadata) for seq in jsonpangenome.sequences},
                    )
    seq_int_to_str_ids = {s.sequence_int_id: s.sequence_str_id for s in jsonpangenome.sequences}
    consensus_tree = ConsensusTree()
    consensus_tree.nodes = [ConsensusNode(consensus_id=ConsensusNodeID(c.consensus_node_id),
                            parent_node_id=ConsensusNodeID(c.parent),
                            children_nodes_ids=[ConsensusNodeID(child) for child in c.children],
                            sequences_ids=[SequenceID(seq_int_to_str_ids[seq_int_id]) for seq_int_id in c.sequences_int_ids],
                            mincomp=CompatibilityToPath(c.mincomp)) for c in jsonpangenome.consensuses]

    return poagraph, consensus_tree

def expand_unassigned_node(consensus_tree: ConsensusTree):
    node_to_remove = None
    new_nodes = []
    for node in consensus_tree.nodes:
        if node.parent_node_id == ConsensusNodeID(0) and node.mincomp == CompatibilityToPath(0):
            node_to_remove = node.consensus_id
            for sequence in node.sequences_ids:
                new_nodes.append(ConsensusNode(consensus_id=ConsensusNodeID(len(consensus_tree.nodes)+len(new_nodes)),
                                               parent_node_id=ConsensusNodeID(0),
                                               children_nodes_ids=[],
                                               sequences_ids=[sequence],
                                               mincomp=CompatibilityToPath(1.0),
                                               compatibilities_to_all={},
                                               consensus_path=[]))
            node.sequences_ids = []
            node.children_nodes_ids = []
        elif node.consensus_id != ConsensusNodeID(0):
            node.children_nodes_ids = []
    for new_node in new_nodes:
        consensus_tree.nodes[0].children_nodes_ids.append(new_node.consensus_id)
    consensus_tree.nodes.extend(new_nodes)
    return consensus_tree

for sample in range(0, 1):
    for p in ["08", "099"]:
    # for p in ["0995"]:
    #     dir_path = f"../data/Simulated3_fewcycles/output/analysis_1508_105936/{sample}_{p}__output/"
        dir_path = f"../data/Ebola/output/analysis_1408_034144/{sample}_{p}__output/"
        # dir_path = f"../data/Simulated1/output/analysis_0808_035057/{sample}_{p}__output/"
        poagraph, consensus_tree = convert_jsonpangenome(path=Path(f"{dir_path}/pangenome.json"))
        # seq_id_to_metadata = {seq_id: seq.seqmetadata for seq_id, seq in poagraph.sequences.items()}
        seq_id_to_metadata = None
        consensus_tree = expand_unassigned_node(consensus_tree)
        newick_consensus_tree = consensus_tree.as_newick(seq_id_to_metadata, expand_leaves=True)
        pathtools.save_to_file(newick_consensus_tree, Path(f"{dir_path}/consensus_tree_extended_2.newick"))
