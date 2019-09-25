from collections import deque
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np

from pangtreebuild.affinitytree.AffinityTree import AffinityTree, AffinityNode, AffinityNodeID, Compatibility
from pangtreebuild.affinitytree.parameters import Blosum, Stop, P, Hbmin
from pangtreebuild.datamodel.Poagraph import Poagraph
from pangtreebuild.affinitytree import poa
from pangtreebuild.datamodel.Sequence import SequenceID, SeqPath
from pangtreebuild.tools import logprocess

tresholds_logger = logprocess.get_logger('tresholdsCSV')
detailed_logger = logprocess.get_logger('details')
global_logger = logprocess.get_global_logger()


class AffinityTreeNodeGenerationException(Exception):
    pass


def get_affinity_tree(poagraph: Poagraph,
                      blosum: Blosum,
                      output_dir: Path,
                      stop: Stop,
                      p: P,
                      verbose: bool) -> AffinityTree:
    global_logger.info("Affinity Tree generation started.")
    if verbose:
        logprocess.add_file_handler_to_logger(output_dir, "tresholdsCSV", "tresholds.csv", "%(message)s", False)
    _raise_error_if_invalid_poagraph(poagraph)
    affinity_tree = _init_affinity_tree(poagraph, blosum.filepath, output_dir, p)

    nodes_to_process = deque([affinity_tree.get_node(AffinityNodeID(0))])
    while nodes_to_process:
        node = nodes_to_process.pop()

        children_nodes = _get_children_nodes_looping(node,
                                            poagraph,
                                            output_dir,
                                            blosum.filepath,
                                            p,
                                            affinity_tree.get_max_node_id())
        if len(children_nodes) == 1:
            continue

        for child in children_nodes:
            child.compatibilities = poagraph.get_compatibilities(sequences_ids=[*poagraph.sequences.keys()],
                                                                 consensus_path=child.consensus,
                                                                 p=p)
            node.children.append(child.id_)
            affinity_tree.nodes.append(child)
            if not _node_is_ready(child, stop):
                nodes_to_process.append(child)
    global_logger.info("Affinity Tree generation finished.\n")
    return affinity_tree


def _raise_error_if_invalid_poagraph(poagraph: Poagraph):
    if len(poagraph.sequences) == 0:
        raise AffinityTreeNodeGenerationException("Invalid poagraph."
                                                  "No paths in poagraph."
                                                  "Affinity Tree generation is impossible.")


def _init_affinity_tree(poagraph: Poagraph, blosum_path: Path, output_dir: Path, p: P) -> AffinityTree:
    affinity_tree = AffinityTree()
    root_node = _get_root_node(poagraph, blosum_path, output_dir, p)
    affinity_tree.nodes.append(root_node)
    return affinity_tree


def _get_root_node(poagraph: Poagraph, blosum_path: Path, output_dir: Path, p: P) -> AffinityNode:
    detailed_logger.info("Getting the root affinity node...")
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
        raise AffinityTreeNodeGenerationException("Cannot find root consensus.")
    compatibilities = poagraph.get_compatibilities(all_poagraph_sequences_ids,
                                                   consensus_paths[0].path,
                                                   p=p)
    affinity_node = AffinityNode(id_=AffinityNodeID(0),
                                 sequences=[*poagraph.sequences.keys()],
                                 mincomp=_get_min_comp(all_poagraph_sequences_ids, compatibilities),
                                 compatibilities=compatibilities,
                                 consensus=consensus_paths[0].path)
    detailed_logger.info(f"New affinity node created: {str(affinity_node)}")
    return affinity_node


def _get_min_comp(node_sequences_ids: List[SequenceID],
                  comps_to_consensus: Dict[SequenceID, Compatibility]) -> Compatibility:
    compatibilities_of_node_sequences = [comp
                                         for seq_id, comp
                                         in comps_to_consensus.items()
                                         if seq_id in node_sequences_ids]
    if not compatibilities_of_node_sequences:
        raise AffinityTreeNodeGenerationException("Cannot provide mincomp."
                                               "No sequences assigned to this affinity node.")
    return min(compatibilities_of_node_sequences)


def _get_children_nodes_looping(node: AffinityNode,
                                poagraph: Poagraph,
                                output_dir: Path,
                                blosum_path: Path,
                                p: P,
                                current_max_affinity_node_id: int) -> List[AffinityNode]:
    children_nodes: List[AffinityNode] = []
    not_assigned_sequences_ids: List[SequenceID] = node.sequences
    detailed_logger.info(f"Getting children nodes for affinity node {node.id_}...")

    affinity_node_id = 0
    so_far_cutoffs: List[Compatibility] = []
    while not_assigned_sequences_ids:
        detailed_logger.info(f"### Getting child {len(so_far_cutoffs)}...")
        child_ready = False
        attempt = 0
        current_candidates = not_assigned_sequences_ids
        while not child_ready:
            consensus_candidate = poa.get_consensuses(poagraph,
                                                      current_candidates,
                                                      output_dir,
                                                      f"parent_{node.id_}_child_{len(so_far_cutoffs)}_attempt_{attempt}",
                                                      blosum_path,
                                                      Hbmin(0),
                                                      specific_consensuses_id=[0])[0].path
            compatibilities_to_consensus_candidate = poagraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                                  consensus_path=consensus_candidate,
                                                                                  p=p)
            compatibilities_to_consensus_candidate[SequenceID("parent")] = node.mincomp
            qualified_sequences_ids_candidates, cutoff = _get_qualified_sequences_ids_and_cutoff(
                compatibilities_to_max_c=compatibilities_to_consensus_candidate,
                so_far_cutoffs=so_far_cutoffs,
                splitted_node_id=node.id_)

            if qualified_sequences_ids_candidates == current_candidates or attempt == 10:
                if attempt == 10:
                    detailed_logger.info("Attempt treshold 10 exceeded!")
                affinity_node_id += 1

                affinity_node = AffinityNode(
                    id_=AffinityNodeID(current_max_affinity_node_id + affinity_node_id),
                    parent=node.id_,
                    sequences=qualified_sequences_ids_candidates,
                    mincomp=_get_min_comp(node_sequences_ids=qualified_sequences_ids_candidates,
                                          comps_to_consensus=compatibilities_to_consensus_candidate),
                    consensus=SeqPath(consensus_candidate))
                children_nodes.append(affinity_node)
                not_assigned_sequences_ids = list(set(not_assigned_sequences_ids) - set(qualified_sequences_ids_candidates))
                child_ready = True
                so_far_cutoffs.append(affinity_node.mincomp)
            else:
                current_candidates = qualified_sequences_ids_candidates
                attempt += 1

    detailed_logger.info("Children nodes generated.")

    return children_nodes

#
# def _get_children_nodes(node: ConsensusNode,
#                         poagraph: Poagraph,
#                         output_dir: Path,
#                         blosum_path: Path,
#                         p: P,
#                         max_cutoff_strategy: FindMaxCutoff,
#                         node_cutoff_strategy: FindNodeCutoff,
#                         current_max_affinity_node_id: int) -> List[ConsensusNode]:
#         children_nodes: List[ConsensusNode] = []
#         not_assigned_sequences_ids: List[SequenceID] = node.sequences
#         so_far_cutoffs: List[CompatibilityToPath] = []
#         detailed_logger.info(f"Getting children nodes for consensus node {node.id_}...")
#
#         while not_assigned_sequences_ids:
#             detailed_logger.info(f"### Getting child {len(so_far_cutoffs)}...")
#             consensus = poa.get_consensuses(poagraph,
#                                                  not_assigned_sequences_ids,
#                                                  output_dir,
#                                                  f"{node.id_}_{len(so_far_cutoffs)}_all",
#                                                  blosum_path,
#                                                  Hbmin(0),
#                                                  specific_consensuses_id=[0])[0].path
#
#             compatibilities_to_consensus = poagraph.get_compatibilities(sequences=not_assigned_sequences_ids,
#                                                                         consensus=consensus,
#                                                                         p=p)
#             compatibilities_to_consensus[SequenceID("parent")] = node.mincomp
#
#             max_sequences_ids, _ = _get_max_compatible_sequences_ids_and_cutoff(max_cutoff_strategy,
#                                                                   compatibilities_to_consensus,
#                                                                   splitted_node_id=node.id_)
#             max_consensus_path = poa.get_consensuses(poagraph,
#                                                      max_sequences_ids,
#                                                      output_dir,
#                                                      f"{node.id_}_{len(so_far_cutoffs)}_max",
#                                                      blosum_path,
#                                                      Hbmin(0),
#                                                      specific_consensuses_id=[0])[0].path
#
#             comps_to_max_consensus = poagraph.get_compatibilities(sequences=not_assigned_sequences_ids,
#                                                                   consensus=max_consensus_path,
#                                                                   p=p)
#             comps_to_max_consensus[SequenceID("parent")] = node.mincomp
#
#             qualified_sequences_ids, node_cutoff = _get_qualified_sequences_ids_and_cutoff(
#                 node_cutoff_strategy,
#                 compatibilities_to_max_c=comps_to_max_consensus,
#                 so_far_cutoffs=so_far_cutoffs, splitted_node_id=node.id_)
#             consensus_node = ConsensusNode(id_=AffinityNodeID(current_max_affinity_node_id + len(so_far_cutoffs)+1),
#                                            parent=node.id_,
#                                            sequences=qualified_sequences_ids,
#                                            mincomp=_get_min_comp(node_sequences_ids=qualified_sequences_ids,
#                                                                  comps_to_consensus=comps_to_max_consensus),
#                                            consensus=Path(max_consensus_path))
#             detailed_logger.info(f"New consensus node created: {str(consensus_node)}")
#             so_far_cutoffs.append(node_cutoff)
#             children_nodes.append(consensus_node)
#             not_assigned_sequences_ids = list(set(not_assigned_sequences_ids) - set(qualified_sequences_ids))
#
#         detailed_logger.info("Children nodes generated.")
#
#         return children_nodes

#
# def _get_max_compatible_sequences_ids_and_cutoff(compatibilities_to_consensus: Dict[SequenceID, CompatibilityToPath],
#                                                  splitted_node_id,
#                                                  so_far_cutoffs: List[CompatibilityToPath] = []) -> List[SequenceID]:
#     max_cutoff = max_cutoff_strategy.find_max_distance([*compatibilities_to_consensus.values()])
#
#     tresholds_logger.info(f"Splitting {splitted_node_id}; MAX; {compatibilities_to_consensus}; "
#                           f"{max_cutoff.cutoff}; {max_cutoff.explanation}")
#     detailed_logger.info(f"Minimum compatibility of sequences chosen for generating the best consensus: "
#                          f"{max_cutoff.cutoff}.")
#     max_sequences_ids = _get_sequences_ids_above_cutoff(compatibilities_to_consensus, max_cutoff.cutoff)
#     return max_sequences_ids, max_cutoff


def _get_sequences_ids_above_cutoff(compatibilities: Dict[SequenceID, Compatibility],
                                    cutoff: Compatibility) -> List[SequenceID]:
    return [seq_id for seq_id, comp in compatibilities.items() if comp >= cutoff and seq_id != SequenceID("parent")]


def _get_qualified_sequences_ids_and_cutoff(compatibilities_to_max_c: Dict[SequenceID, Compatibility],
                                            so_far_cutoffs: List[Compatibility], splitted_node_id) \
        -> Tuple[List[SequenceID], Compatibility]:
    node_cutoff = find_node_cutoff(compatibilities=[*compatibilities_to_max_c.values()],
                                                        so_far_cutoffs=so_far_cutoffs)
    tresholds_logger.info(f"Splitting {splitted_node_id}; NODE; {compatibilities_to_max_c}; "
                          f"{node_cutoff.cutoff}; {node_cutoff.explanation}")
    compatibtle_sequences_ids = _get_sequences_ids_above_cutoff(compatibilities_to_max_c, node_cutoff.cutoff)
    detailed_logger.info(f"{len(compatibtle_sequences_ids)} sequences ({compatibtle_sequences_ids} are "
                         f"qualified to be enclosed in this node, mincomp = {node_cutoff.cutoff}.")
    return compatibtle_sequences_ids, node_cutoff.cutoff


def _node_is_ready(node: AffinityNode, stop: Stop) -> bool:
    if len(node.sequences) == 1 or node.mincomp.base_value() >= stop:
        detailed_logger.info(f"Node {node.id_} satisfied requirements and won't be split!")
        return True
    return False


class FindCutoffResult:
    def __init__(self, cutoff: Compatibility, explanation: str):
        self.cutoff: Compatibility = cutoff
        self.explanation: str = explanation

def find_node_cutoff(compatibilities: List[Compatibility],
                     so_far_cutoffs: List[Compatibility]) -> FindCutoffResult:
    if not so_far_cutoffs:
        cutoff = find_max_distance(compatibilities)
        reason = "No so far cutoffs. Find max distance in sorted values."
    else:
        guard = min(so_far_cutoffs)
        sorted_comp = sorted(compatibilities)
        if guard <= sorted_comp[0]:
            cutoff = sorted_comp[0]
            reason = "guard < min(compatibilities). Return min(compatibilities)."
        elif guard >= sorted_comp[-1]:
            cutoff = find_max_distance(compatibilities)
            reason = "guard > max(compatibilities). Find max distance in sorted values."
        else:
            first_comp_greater_than_guard_index = [i for i, c in enumerate(sorted_comp) if c > guard][0]
            cutoff = find_max_distance(sorted_comp[0:first_comp_greater_than_guard_index + 1])
            reason = "Find max distance in sorted_comp[0:first_comp_greater_than_guard_index + 1]"
    return FindCutoffResult(cutoff, reason)


def find_max_distance(compatibilities: List[Compatibility]) -> Compatibility:
    def break_if_empty(compatibilities: List[Compatibility]) -> None:
        if not list(compatibilities):
            raise ValueError("Empty compatibilities list. Cannot find cutoff.")

    break_if_empty(compatibilities)
    if len(compatibilities) == 1:
        return compatibilities[0]

    sorted_values = sorted(compatibilities)
    distances = np.array([sorted_values[i + 1] - sorted_values[i] for i in range(len(sorted_values) - 1)])
    max_distance_index: int = np.argmax(distances)
    return sorted_values[max_distance_index + 1]





