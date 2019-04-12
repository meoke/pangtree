from pathlib import Path
from typing import Dict

from consensus import poa
from consensus.ConsensusTree import ConsensusTree, ConsensusNode, ConsensusNodeID
from consensus.input_types import Blosum, Hbmin
from consensus.poa import ConsInfo
from datamodel.Poagraph import Poagraph
from datamodel.Sequence import SequencePath
from tools import logprocess

global_logger = logprocess.get_global_logger()

class ConsensusGenerationException(Exception):
    pass


def _get_consensus_nodes(poagraph: Poagraph, consensus_paths: Dict[int, ConsInfo]):
    consensus_nodes = []
    assigned_sequences = []
    for c_id, c_info in consensus_paths.items():
        assigned_sequences += c_info.assigned_sequences_ids
        compatibilities = poagraph.get_compatibilities(poagraph.get_sequences_ids(), c_info.path)
        mincomp = min([c for seq_id, c in compatibilities.items() if seq_id in c_info.assigned_sequences_ids])
        cn = ConsensusNode(consensus_id=ConsensusNodeID(c_id+1),
                           parent_node_id=ConsensusNodeID(0),
                           sequences_ids=c_info.assigned_sequences_ids,
                           mincomp=mincomp,
                           compatibilities_to_all=compatibilities,
                           consensus_path=c_info.path)
        consensus_nodes.append(cn)

    node_for_unassigned_sequences = ConsensusNode(parent_node_id=ConsensusNodeID(0),
                                                  sequences_ids=[seq_id
                                                                 for seq_id in poagraph.get_sequences_ids()
                                                                 if seq_id not in assigned_sequences],
                                                  consensus_id=ConsensusNodeID(len(consensus_nodes)+1))
    consensus_nodes.append(node_for_unassigned_sequences)
    print(consensus_nodes)
    return consensus_nodes


def get_simple_consensus_tree(poagraph: Poagraph, blosum: Blosum, output_dir: Path, hbmin: Hbmin, verbose: bool):
    global_logger.info("Simple consensus generation started.")
    all_poagraph_sequences_ids = poagraph.get_sequences_ids()
    try:
        consensus_paths = poa.get_consensuses(poagraph,
                                              poagraph.get_sequences_ids(),
                                              output_dir,
                                              "simple_tree",
                                              blosum.filepath,
                                              hbmin)
    except poa.NoConsensusError:
        raise ConsensusGenerationException("Cannot find root consensus.")
    consensus_tree = ConsensusTree()
    consensus_nodes = _get_consensus_nodes(poagraph, consensus_paths)
    root_node = ConsensusNode(consensus_id=ConsensusNodeID(0),
                              children_nodes_ids=[consensus_node.consensus_id
                                                  for consensus_node in consensus_nodes])
    consensus_tree.nodes.append(root_node)
    consensus_tree.nodes += consensus_nodes
    return consensus_tree
