from pathlib import Path
from collections import deque
from typing import List, Dict, Tuple

from consensus.ConsensusesTree import ConsensusesTree
from consensus.ConsensusNode import ConsensusNode, Compatibility, ConsensusNodeID
from consensus.FindCutoff import FindMaxCutoff, FindNodeCutoff
from consensus.exceptions import TreeConsensusGenerationException
from pangraph.Pangraph import Pangraph
from pangraph.custom_types import SequenceID, NodeID
import consensus.top_consensus as top_consensus
from tools import loggingtools

tresholds_logger = loggingtools.get_logger('tresholdsCSV')
detailed_logger = loggingtools.get_logger('details')
global_logger = loggingtools.get_global_logger()

class TreePOAConsensusGenerator:
    def __init__(self,
                 max_node_strategy: FindMaxCutoff,
                 node_cutoff_strategy: FindNodeCutoff,
                 blosum_path: Path,
                 stop: Compatibility,
                 re_consensus: bool):
        self.max_node_strategy: FindMaxCutoff = max_node_strategy
        self.node_cutoff_strategy: FindNodeCutoff = node_cutoff_strategy
        self.blosum_path: Path = blosum_path
        self.stop: Compatibility = stop
        self.re_consensus: bool = re_consensus
        self.pangraph: Pangraph = None
        self.output_dir: Path = None
        self.consensuses_tree: ConsensusesTree = None



    def get_consensuses_tree(self,
                             pangraph: Pangraph,
                             output_dir: Path,
                             log_tresholds: bool
                             ) -> ConsensusesTree:
        global_logger.info("Consensuses Tree generation started.")
        self.pangraph = pangraph
        self.output_dir = output_dir
        # self.set_tresholds_log_handler()
        if log_tresholds:
            loggingtools.add_fileHandler_to_logger(output_dir, "tresholdsCSV", "tresholds.csv", "%(message)s", False)

        self.break_if_pangraph_is_invalid()

        self.consensuses_tree = self.init_consensuses_tree()

        nodes_to_process = deque([self.consensuses_tree.get_node(ConsensusNodeID(0))])
        while nodes_to_process:
            node = nodes_to_process.pop()
            children_nodes = self.get_children_nodes(node)

            if len(children_nodes) == 1:
                continue

            for child in children_nodes:
                child.compatibilities_to_all = pangraph.get_compatibilities(sequences_ids=[*pangraph.paths.keys()],
                                                                            consensus=child.consensus_path)
                node.children_nodes_ids.append(child.consensus_id)
                self.consensuses_tree.nodes.append(child)
                if not self.node_is_ready(child):
                    nodes_to_process.append(child)
        global_logger.info("Consensuses Tree generation finished.\n")
        return self.consensuses_tree

    def init_consensuses_tree(self) -> ConsensusesTree:
        consensuses_tree = ConsensusesTree()
        root_node = self.get_root_node()
        consensuses_tree.nodes.append(root_node)
        return consensuses_tree

    def get_top_consensus(self, sequences_ids: List[SequenceID], name: str) -> List[NodeID]:
        try:
            return top_consensus.get_top_consensus(pangraph=self.pangraph,
                                                   sequences_ids=sequences_ids,
                                                   output_dir=self.output_dir,
                                                   file_prefix=name,
                                                   blosum_path=self.blosum_path)
        except TreeConsensusGenerationException as e:
            raise TreeConsensusGenerationException(f'Cannot find {name} consensus.') from e

    def get_root_node(self) -> ConsensusNode:
        detailed_logger.info("Getting the root consensus node...")
        consensus_path = self.get_top_consensus(sequences_ids=[*self.pangraph.paths.keys()],
                                                name="root")
        all_pangraph_sequences_ids = [*self.pangraph.paths.keys()]
        compatibilities = self.pangraph.get_compatibilities(sequences_ids=all_pangraph_sequences_ids,
                                                            consensus=consensus_path)
        consensus_node = ConsensusNode(sequences_ids=all_pangraph_sequences_ids,
                                       consensus_id=ConsensusNodeID(0),
                                       mincomp=self.get_mincomp(all_pangraph_sequences_ids, compatibilities),
                                       compatibilities_to_all=compatibilities,
                                       consensus_path=consensus_path)
        detailed_logger.info(f"New consensus node created: {str(consensus_node)}")
        return consensus_node

    def break_if_pangraph_is_invalid(self) -> None:
        if len(self.pangraph.paths) == 0:
            raise TreeConsensusGenerationException("Invalid pangraph."
                                                   "No paths in pangraph."
                                                   "Cannot find consensuses.")

    def node_is_ready(self, node: ConsensusNode) -> bool:
        if len(node.sequences_ids) == 1 or node.mincomp >= self.stop:
            detailed_logger.info(f"Node {node.consensus_id} satisfied requirements and won't be split!")
            return True
        return False

    def get_children_nodes(self, node: ConsensusNode) -> List[ConsensusNode]:
        children_nodes: List[ConsensusNode] = []
        not_assigned_sequences_ids: List[SequenceID] = node.sequences_ids
        so_far_cutoffs: List[Compatibility] = []
        detailed_logger.info(f"Getting children nodes for consensus node {node.consensus_id}...")

        while not_assigned_sequences_ids:
            detailed_logger.info(f"### Getting child {len(so_far_cutoffs)}...")
            consensus_path = self.get_top_consensus(sequences_ids=not_assigned_sequences_ids,
                                                    name=f"{node.consensus_id}_{len(so_far_cutoffs)}_all")
            compatibilities_to_consensus = self.pangraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                             consensus=consensus_path)

            max_sequences_ids = self.get_max_compatible_sequences_ids(compatibilities_to_consensus, splitted_node_id=node.consensus_id)
            max_consensus_path = self.get_top_consensus(sequences_ids=max_sequences_ids,
                                                        name=f"{node.consensus_id}_{len(so_far_cutoffs)}_max")

            comps_to_max_consensus = self.pangraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                       consensus=max_consensus_path)
            qualified_sequences_ids, node_cutoff = self.get_qualified_sequences_ids_and_cutoff(
                compatibilities_to_max_c=comps_to_max_consensus,
                so_far_cutoffs=so_far_cutoffs, splitted_node_id=node.consensus_id)
            consensus_node = ConsensusNode(parent_node_id=node.consensus_id,
                                           sequences_ids=qualified_sequences_ids,
                                           consensus_id=ConsensusNodeID(len(self.consensuses_tree.nodes) + len(so_far_cutoffs)),
                                           mincomp=self.get_mincomp(node_sequences_ids=qualified_sequences_ids,
                                                                    comps_to_consensus=comps_to_max_consensus),
                                           consensus_path=max_consensus_path)
            detailed_logger.info(f"New consensus node created: {str(consensus_node)}")
            so_far_cutoffs.append(node_cutoff)
            children_nodes.append(consensus_node)
            not_assigned_sequences_ids = list(set(not_assigned_sequences_ids) - set(qualified_sequences_ids))

        if self.re_consensus:
            detailed_logger.info("Trying to re consensus...")
            children_nodes = self.reorder_consensuses(children_nodes)
        detailed_logger.info("Children nodes generated.")

        return children_nodes

    def get_max_compatible_sequences_ids(self, compatibilities_to_consensus: Dict[SequenceID, Compatibility], splitted_node_id) ->\
            List[SequenceID]:
        max_cutoff = self.max_node_strategy.find_max_cutoff([*compatibilities_to_consensus.values()])
        tresholds_logger.info(f"Splitting {splitted_node_id}; MAX; {compatibilities_to_consensus}; {max_cutoff.cutoff}; {max_cutoff.explanation}")
        detailed_logger.info(f"Minimum compatibility of sequences chosen for generating the best consensus: {max_cutoff.cutoff}.")
        max_sequences_ids = self.get_sequences_ids_above_cutoff(compatibilities_to_consensus, max_cutoff.cutoff)
        return max_sequences_ids

    def get_qualified_sequences_ids_and_cutoff(self,
                                               compatibilities_to_max_c: Dict[SequenceID, Compatibility],
                                               so_far_cutoffs: List[Compatibility], splitted_node_id) \
            -> Tuple[List[SequenceID], Compatibility]:
        node_cutoff = self.node_cutoff_strategy.find_node_cutoff(compatibilities=[*compatibilities_to_max_c.values()],
                                                                 so_far_cutoffs=so_far_cutoffs)
        tresholds_logger.info(f"Splitting {splitted_node_id}; NODE; {compatibilities_to_max_c}; {node_cutoff.cutoff}; {node_cutoff.explanation}")
        compatibtle_sequences_ids = self.get_sequences_ids_above_cutoff(compatibilities_to_max_c, node_cutoff.cutoff)
        detailed_logger.info(f"{len(compatibtle_sequences_ids)} sequences ({compatibtle_sequences_ids} are "
                     f"qualified to be enclosed in this node, mincomp = {node_cutoff.cutoff}.")
        return compatibtle_sequences_ids, node_cutoff.cutoff

    def get_mincomp(self,
                    node_sequences_ids: List[SequenceID],
                    comps_to_consensus: Dict[SequenceID, Compatibility]) -> Compatibility:
        compatibilities_of_node_sequences = [comp
                                             for seq_id, comp
                                             in comps_to_consensus.items()
                                             if seq_id in node_sequences_ids]
        if not compatibilities_of_node_sequences:
            raise TreeConsensusGenerationException("Cannot provide mincomp."
                                                   "No sequences assigned to this consensus node.")
        return min(compatibilities_of_node_sequences)

    def get_sequences_ids_above_cutoff(self,
                                       compatibilities: Dict[SequenceID, Compatibility],
                                       cutoff: Compatibility) ->\
            List[SequenceID]:
        return [seq_id for seq_id, comp in compatibilities.items() if comp >= cutoff]

    def reorder_consensuses(self, children_nodes: List[ConsensusNode]) -> List[ConsensusNode]:
        any_sequence_was_moved = False
        all_sequences_ids = [seq_id for c_node in children_nodes for seq_id in c_node.sequences_ids]
        possible_compatibilities = {seq_id: [] for seq_id in all_sequences_ids}
        current_compatibilities_assignment = {}
        for consensus in children_nodes:
            seq_id_to_comp = self.pangraph.get_compatibilities(all_sequences_ids, consensus.consensus_path)
            for seq_id, comp in seq_id_to_comp.items():
                possible_compatibilities[seq_id].append((consensus.consensus_id, comp))
            for seq_id in consensus.sequences_ids:
                current_compatibilities_assignment[seq_id] = (consensus.consensus_id, seq_id_to_comp[seq_id])

        for seq_id, comp_id_comp in possible_compatibilities.items():
            sorted_possible_compatibilities = sorted(comp_id_comp, key=lambda comp_id_comp: comp_id_comp[1], reverse=True)
            current_comp = current_compatibilities_assignment[seq_id][1]
            current_comp_id = current_compatibilities_assignment[seq_id][0]
            if sorted_possible_compatibilities[0][1] > current_comp:
                better_comp_id = sorted_possible_compatibilities[0][0]
                children_nodes = self.move(children_nodes, current_comp_id=current_comp_id, better_comp_id=better_comp_id, seq_id=seq_id)
                any_sequence_was_moved = True
        if not any_sequence_was_moved:
            detailed_logger.info("No sequence was reassigned to different node.")

        return self.update_children_nodes(children_nodes)

    def move(self, children_nodes, current_comp_id, better_comp_id, seq_id):
        detailed_logger.info(f"Sequence {seq_id} moved from consensus {current_comp_id} to {better_comp_id}.")
        for consensus in children_nodes:
            if consensus.consensus_id == current_comp_id:
                to_del = consensus.sequences_ids.index(seq_id)
                del consensus.sequences_ids[to_del]
            elif consensus.consensus_id == better_comp_id:
                consensus.sequences_ids.append(seq_id)
        return children_nodes

    def update_children_nodes(self, children_nodes: List[ConsensusNode]):
        to_remove = []
        for i, child in enumerate(children_nodes):
            if not child.sequences_ids:
                to_remove.append(i)
            else:
                new_compatibilities = self.pangraph.get_compatibilities(child.sequences_ids, child.consensus_path)
                child.mincomp = min([comp for seq_id, comp in new_compatibilities.items()])
        for r in sorted(to_remove, reverse=True):
            del children_nodes[r]
        return children_nodes

