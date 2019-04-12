from collections import deque
from pathlib import Path
from typing import List, Dict, Tuple

from consensus.ConsensusTree import ConsensusTree, ConsensusNode, ConsensusNodeID, CompatibilityToPath
from consensus.cutoffs import FindMaxCutoff, FindNodeCutoff
from consensus.input_types import Blosum, Stop, P, Hbmin
from datamodel.Poagraph import Poagraph
from consensus import poa
from datamodel.Sequence import SequenceID, SequencePath
from tools import logprocess

tresholds_logger = logprocess.get_logger('tresholdsCSV')
detailed_logger = logprocess.get_logger('details')
global_logger = logprocess.get_global_logger()


class TreeConsensusGenerationException(Exception):
    pass


def get_consensus_tree(poagraph: Poagraph,
                       blosum: Blosum,
                       output_dir: Path,
                       stop: Stop,
                       p: P,
                       max_strategy: FindMaxCutoff,
                       node_strategy: FindNodeCutoff,
                       verbose: bool) -> ConsensusTree:
    global_logger.info("Consensuses Tree generation started.")
    if verbose:
        logprocess.add_file_handler_to_logger(output_dir, "tresholdsCSV", "tresholds.csv", "%(message)s", False)
    _raise_error_if_invalid_poagraph(poagraph)
    consensus_tree = _init_consensus_tree(poagraph, blosum.filepath, output_dir, p)

    nodes_to_process = deque([consensus_tree.get_node(ConsensusNodeID(0))])
    while nodes_to_process:
        node = nodes_to_process.pop()

        children_nodes = _get_children_nodes(node,
                                            poagraph,
                                            output_dir,
                                            blosum.filepath,
                                            p,
                                            max_strategy,
                                            node_strategy,
                                            consensus_tree.get_max_node_id())
        if len(children_nodes) == 1:
            continue

        for child in children_nodes:
            child.compatibilities_to_all = poagraph.get_compatibilities(sequences_ids=[*poagraph.sequences.keys()],
                                                                        consensus_path=child.consensus_path,
                                                                        p=p)
            node.children_nodes_ids.append(child.consensus_id)
            consensus_tree.nodes.append(child)
            if not _node_is_ready(child, stop):
                nodes_to_process.append(child)
    global_logger.info("Consensuses Tree generation finished.\n")
    return consensus_tree


def _raise_error_if_invalid_poagraph(poagraph: Poagraph):
    if len(poagraph.sequences) == 0:
        raise TreeConsensusGenerationException("Invalid pangraph."
                                               "No paths in pangraph."
                                               "Cannot find consensuses.")


def _init_consensus_tree(poagraph: Poagraph, blosum_path: Path, output_dir: Path, p: P) -> ConsensusTree:
    consensuses_tree = ConsensusTree()
    root_node = _get_root_node(poagraph, blosum_path, output_dir, p)
    consensuses_tree.nodes.append(root_node)
    return consensuses_tree


def _get_root_node(poagraph: Poagraph, blosum_path: Path, output_dir: Path, p: P) -> ConsensusNode:
    detailed_logger.info("Getting the root consensus node...")
    all_poagraph_sequences_ids = poagraph.get_sequences_ids()
    try:
        consensus_paths = poa.get_consensuses(poagraph,
                                              all_poagraph_sequences_ids,
                                              output_dir,
                                              "root",
                                              blosum_path,
                                              hbmin=Hbmin(0),
                                              specific_consensuses_id=[0])
    except poa.NoConsensusError:
        raise TreeConsensusGenerationException("Cannot find root consensus.")
    compatibilities = poagraph.get_compatibilities(all_poagraph_sequences_ids,
                                                   consensus_paths[0].path,
                                                   p=p)
    consensus_node = ConsensusNode(consensus_id=ConsensusNodeID(0),
                                   sequences_ids=[*poagraph.sequences.keys()],
                                   mincomp=_get_min_comp(all_poagraph_sequences_ids, compatibilities),
                                   compatibilities_to_all=compatibilities,
                                   consensus_path=consensus_paths[0].path)
    detailed_logger.info(f"New consensus node created: {str(consensus_node)}")
    return consensus_node


def _get_min_comp(node_sequences_ids: List[SequenceID],
                  comps_to_consensus: Dict[SequenceID, CompatibilityToPath]) -> CompatibilityToPath:
    compatibilities_of_node_sequences = [comp
                                         for seq_id, comp
                                         in comps_to_consensus.items()
                                         if seq_id in node_sequences_ids]
    if not compatibilities_of_node_sequences:
        raise TreeConsensusGenerationException("Cannot provide mincomp."
                                               "No sequences assigned to this consensus node.")
    return min(compatibilities_of_node_sequences)


def _get_children_nodes(node: ConsensusNode,
                        poagraph: Poagraph,
                        output_dir: Path,
                        blosum_path: Path,
                        p: P,
                        max_cutoff_strategy: FindMaxCutoff,
                        node_cutoff_strategy: FindNodeCutoff,
                        current_max_consensus_node_id: int) -> List[ConsensusNode]:
        children_nodes: List[ConsensusNode] = []
        not_assigned_sequences_ids: List[SequenceID] = node.sequences_ids
        so_far_cutoffs: List[CompatibilityToPath] = []
        detailed_logger.info(f"Getting children nodes for consensus node {node.consensus_id}...")

        while not_assigned_sequences_ids:
            detailed_logger.info(f"### Getting child {len(so_far_cutoffs)}...")
            consensus_path = poa.get_consensuses(poagraph,
                                                 not_assigned_sequences_ids,
                                                 output_dir,
                                                 f"{node.consensus_id}_{len(so_far_cutoffs)}_all",
                                                 blosum_path,
                                                 Hbmin(0),
                                                 specific_consensuses_id=[0])[0].path

            compatibilities_to_consensus = poagraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                        consensus_path=consensus_path,
                                                                        p=p)

            max_sequences_ids = _get_max_compatible_sequences_ids(max_cutoff_strategy,
                                                                  compatibilities_to_consensus,
                                                                  splitted_node_id=node.consensus_id)
            max_consensus_path = poa.get_consensuses(poagraph,
                                                     max_sequences_ids,
                                                     output_dir,
                                                     f"{node.consensus_id}_{len(so_far_cutoffs)}_max",
                                                     blosum_path,
                                                     Hbmin(0),
                                                     specific_consensuses_id=[0])[0].path

            comps_to_max_consensus = poagraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                  consensus_path=max_consensus_path,
                                                                  p=p)
            qualified_sequences_ids, node_cutoff = _get_qualified_sequences_ids_and_cutoff(
                node_cutoff_strategy,
                compatibilities_to_max_c=comps_to_max_consensus,
                so_far_cutoffs=so_far_cutoffs, splitted_node_id=node.consensus_id)
            consensus_node = ConsensusNode(consensus_id=ConsensusNodeID(current_max_consensus_node_id + len(so_far_cutoffs)+1),
                                           parent_node_id=node.consensus_id,
                                           sequences_ids=qualified_sequences_ids,
                                           mincomp=_get_min_comp(node_sequences_ids=qualified_sequences_ids,
                                                                 comps_to_consensus=comps_to_max_consensus),
                                           consensus_path=SequencePath(max_consensus_path))
            detailed_logger.info(f"New consensus node created: {str(consensus_node)}")
            so_far_cutoffs.append(node_cutoff)
            children_nodes.append(consensus_node)
            not_assigned_sequences_ids = list(set(not_assigned_sequences_ids) - set(qualified_sequences_ids))

        detailed_logger.info("Children nodes generated.")

        return children_nodes


def _get_max_compatible_sequences_ids(max_cutoff_strategy: FindMaxCutoff,
                                      compatibilities_to_consensus: Dict[SequenceID, CompatibilityToPath],
                                      splitted_node_id) -> List[SequenceID]:
    max_cutoff = max_cutoff_strategy.find_max_cutoff([*compatibilities_to_consensus.values()])
    tresholds_logger.info(f"Splitting {splitted_node_id}; MAX; {compatibilities_to_consensus}; "
                          f"{max_cutoff.cutoff}; {max_cutoff.explanation}")
    detailed_logger.info(f"Minimum compatibility of sequences chosen for generating the best consensus: "
                         f"{max_cutoff.cutoff}.")
    max_sequences_ids = _get_sequences_ids_above_cutoff(compatibilities_to_consensus, max_cutoff.cutoff)
    return max_sequences_ids


def _get_sequences_ids_above_cutoff(compatibilities: Dict[SequenceID, CompatibilityToPath],
                                    cutoff: CompatibilityToPath) -> List[SequenceID]:
    return [seq_id for seq_id, comp in compatibilities.items() if comp >= cutoff]


def _get_qualified_sequences_ids_and_cutoff(node_cutoff_strategy: FindNodeCutoff,
                                           compatibilities_to_max_c: Dict[SequenceID, CompatibilityToPath],
                                           so_far_cutoffs: List[CompatibilityToPath], splitted_node_id) \
        -> Tuple[List[SequenceID], CompatibilityToPath]:
    node_cutoff = node_cutoff_strategy.find_node_cutoff(compatibilities=[*compatibilities_to_max_c.values()],
                                                        so_far_cutoffs=so_far_cutoffs)
    tresholds_logger.info(f"Splitting {splitted_node_id}; NODE; {compatibilities_to_max_c}; "
                          f"{node_cutoff.cutoff}; {node_cutoff.explanation}")
    compatibtle_sequences_ids = _get_sequences_ids_above_cutoff(compatibilities_to_max_c, node_cutoff.cutoff)
    detailed_logger.info(f"{len(compatibtle_sequences_ids)} sequences ({compatibtle_sequences_ids} are "
                         f"qualified to be enclosed in this node, mincomp = {node_cutoff.cutoff}.")
    return compatibtle_sequences_ids, node_cutoff.cutoff


def _node_is_ready(node: ConsensusNode, stop: Stop) -> bool:
    if len(node.sequences_ids) == 1 or node.mincomp >= stop:
        detailed_logger.info(f"Node {node.consensus_id} satisfied requirements and won't be split!")
        return True
    return False