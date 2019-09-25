from pathlib import Path
from typing import Dict

from pangtreebuild.affinitytree import poa
from pangtreebuild.affinitytree.AffinityTree import AffinityTree, AffinityNode, AffinityNodeID, Compatibility
from pangtreebuild.affinitytree.parameters import Blosum, Hbmin
from pangtreebuild.affinitytree.poa import ConsInfo
from pangtreebuild.datamodel.Poagraph import Poagraph
from pangtreebuild.tools import logprocess

global_logger = logprocess.get_global_logger()


class ConsensusGenerationException(Exception):
    pass


def _get_poa_affinity_tree_nodes(poagraph: Poagraph, consensus_paths: Dict[int, ConsInfo]):
    consensus_nodes = []
    assigned_sequences = []
    for c_id, c_info in consensus_paths.items():
        assigned_sequences += c_info.assigned_sequences_ids
        compatibilities = poagraph.get_compatibilities(poagraph.get_sequences_ids(), c_info.path)
        if len(c_info.assigned_sequences_ids):
            mincomp = min([c for seq_id, c in compatibilities.items() if seq_id in c_info.assigned_sequences_ids])
        else:
            mincomp = 0
        cn = AffinityNode(id_=AffinityNodeID(c_id + 1),
                          parent=AffinityNodeID(0),
                          sequences=c_info.assigned_sequences_ids,
                          mincomp=mincomp,
                          compatibilities=compatibilities,
                          consensus=c_info.path,
                          children=[])
        consensus_nodes.append(cn)

    node_for_unassigned_sequences = AffinityNode(parent=AffinityNodeID(0),
                                                 sequences=[seq_id
                                                            for seq_id in poagraph.get_sequences_ids()
                                                            if seq_id not in assigned_sequences],
                                                 id_=AffinityNodeID(len(consensus_nodes) + 1),
                                                 mincomp=Compatibility(0),
                                                 children=[])
    consensus_nodes.append(node_for_unassigned_sequences)
    return consensus_nodes


def get_simple_affinity_tree(poagraph: Poagraph, blosum: Blosum, output_dir: Path, hbmin: Hbmin, verbose: bool):
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
    affinity_tree = AffinityTree()
    consensus_nodes = _get_poa_affinity_tree_nodes(poagraph, consensus_paths)
    root_node = AffinityNode(id_=AffinityNodeID(0),
                             children=[consensus_node.id_
                                       for consensus_node in consensus_nodes])
    affinity_tree.nodes.append(root_node)
    affinity_tree.nodes += consensus_nodes
    return affinity_tree
