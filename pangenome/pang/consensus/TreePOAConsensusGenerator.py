from pathlib import Path
from collections import deque
from typing import List, Dict, Tuple
import logging

from consensus.ConsensusesTree import ConsensusesTree
from consensus.ConsensusNode import ConsensusNode, Compatibility, ConsensusNodeID
from consensus.FindCutoff import FindMaxCutoff, FindNodeCutoff
from consensus.exceptions import TreeConsensusGenerationException
from pangraph.Pangraph import Pangraph
from pangraph.custom_types import SequenceID, NodeID
import consensus.top_consensus as top_consensus
from tools import pathtools

logger = logging.getLogger('tresholdsCSV')

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
                             output_dir: Path
                             ) -> ConsensusesTree:
        self.pangraph = pangraph
        self.output_dir = output_dir
        logger.addHandler(self.get_tresholds_log_handler())

        self.break_if_pangraph_is_invalid()

        self.consensuses_tree = self.init_consensuses_tree()

        nodes_to_process = deque([self.consensuses_tree.get_node(ConsensusNodeID(0))])
        while nodes_to_process:
            node = nodes_to_process.pop()
            logger.info("Cos tu sobie napisze,albo nie,albo tak")
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
        consensus_path = self.get_top_consensus(sequences_ids=[*self.pangraph.paths.keys()],
                                                name="root")
        all_pangraph_sequences_ids = [*self.pangraph.paths.keys()]
        compatibilities = self.pangraph.get_compatibilities(sequences_ids=all_pangraph_sequences_ids,
                                                            consensus=consensus_path)
        return ConsensusNode(parent_node_id=None,
                             sequences_ids=all_pangraph_sequences_ids,
                             consensus_id=ConsensusNodeID(0),
                             mincomp=self.get_mincomp(all_pangraph_sequences_ids, compatibilities),
                             compatibilities_to_all=compatibilities,
                             consensus_path=consensus_path)

    def break_if_pangraph_is_invalid(self) -> None:
        if len(self.pangraph.paths) == 0:
            raise TreeConsensusGenerationException("Invalid pangraph."
                                                   "No paths in pangraph."
                                                   "Cannot find consensuses.")

    def node_is_ready(self, node: ConsensusNode) -> bool:
        if len(node.sequences_ids) == 1 or node.mincomp >= self.stop:
            return True
        return False

    def get_children_nodes(self, node: ConsensusNode) -> List[ConsensusNode]:
        children_nodes: List[ConsensusNode] = []
        not_assigned_sequences_ids: List[SequenceID] = node.sequences_ids
        so_far_cutoffs: List[Compatibility] = []

        while not_assigned_sequences_ids:
            consensus_path = self.get_top_consensus(sequences_ids=not_assigned_sequences_ids,
                                                    name=f"{node.consensus_id}_{len(so_far_cutoffs)}_all")
            compatibilities_to_consensus = self.pangraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                             consensus=consensus_path)

            max_sequences_ids = self.get_max_compatible_sequences_ids(compatibilities_to_consensus)
            max_consensus_path = self.get_top_consensus(sequences_ids=max_sequences_ids,
                                                        name=f"{node.consensus_id}_{len(so_far_cutoffs)}_max")
            comps_to_max_consensus = self.pangraph.get_compatibilities(sequences_ids=not_assigned_sequences_ids,
                                                                       consensus=max_consensus_path)

            qualified_sequences_ids, node_cutoff = self.get_qualified_sequences_ids_and_cutoff(
                                                        compatibilities_to_max_c=comps_to_max_consensus,
                                                        so_far_cutoffs=so_far_cutoffs)

            consensus_node = ConsensusNode(parent_node_id=node.consensus_id,
                                           sequences_ids=qualified_sequences_ids,
                                           consensus_id=ConsensusNodeID(len(self.consensuses_tree.nodes) + len(so_far_cutoffs)),
                                           mincomp=self.get_mincomp(node_sequences_ids=qualified_sequences_ids,
                                                                    comps_to_consensus=comps_to_max_consensus),
                                           consensus_path=max_consensus_path)
            so_far_cutoffs.append(node_cutoff)
            children_nodes.append(consensus_node)
            not_assigned_sequences_ids = list(set(not_assigned_sequences_ids) - set(qualified_sequences_ids))

        if self.re_consensus:
            children_nodes = self.reorder_consensuses(children_nodes)

        return children_nodes

    def get_max_compatible_sequences_ids(self, compatibilities_to_consensus: Dict[SequenceID, Compatibility]) ->\
            List[SequenceID]:
        max_cutoff = self.max_node_strategy.find_max_cutoff([*compatibilities_to_consensus.values()], log=True)
        max_sequences_ids = self.get_sequences_ids_above_cutoff(compatibilities_to_consensus, max_cutoff)
        return max_sequences_ids

    def get_qualified_sequences_ids_and_cutoff(self,
                                               compatibilities_to_max_c: Dict[SequenceID, Compatibility],
                                               so_far_cutoffs: List[Compatibility]) \
            -> Tuple[List[SequenceID], Compatibility]:
        node_cutoff = self.node_cutoff_strategy.find_node_cutoff(compatibilities=[*compatibilities_to_max_c.values()],
                                                                 so_far_cutoffs=so_far_cutoffs,
                                                                 log=True)
        compatibtle_sequences_ids = self.get_sequences_ids_above_cutoff(compatibilities_to_max_c, node_cutoff)
        return compatibtle_sequences_ids, node_cutoff

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
        # print("REORDER NOT IMPLEMENTED")
        all_sequences_ids = [seq_id  for c_node in children_nodes for seq_id in c_node.sequences_ids]
        d = {seq_id: [] for seq_id in all_sequences_ids}
        seq_id_to_consensus_id_comp = {}
        for consensus in children_nodes:
            seq_id_to_comp = self.pangraph.get_compatibilities(all_sequences_ids, consensus.consensus_path)
            for seq_id, comp in seq_id_to_comp.items():
                d[seq_id].append((consensus.consensus_id, comp))
            for seq_id in consensus.sequences_ids:
                seq_id_to_consensus_id_comp[seq_id] = (consensus.consensus_id, seq_id_to_comp[seq_id])

        for seq_id, comp_id_comp in d.items():
            sorted_comp_id_comp = sorted(comp_id_comp, key=lambda comp_id_comp: comp_id_comp[1], reverse=True)
            current_comp = seq_id_to_consensus_id_comp[seq_id][1]
            current_comp_id = seq_id_to_consensus_id_comp[seq_id][0]
            if sorted_comp_id_comp[0][1] > current_comp:
                better_comp_id = sorted_comp_id_comp[0][0]
                children_nodes = self.move(children_nodes, current_comp_id=current_comp_id, better_comp_id=better_comp_id, seq_id=seq_id)

        return children_nodes

    def move(self, children_nodes, current_comp_id, better_comp_id, seq_id):
        for consensus in children_nodes:
            if consensus.consensus_id == current_comp_id:
                to_del = consensus.sequences_ids.index(seq_id)
                del consensus.sequences_ids[to_del]
            elif consensus.consensus_id == better_comp_id:
                consensus.sequences_ids.append(seq_id)
        return children_nodes

    def set_tresholds_log_handler(self) -> logging.Handler:
        fh = logging.FileHandler(pathtools.get_child_file_path(self.output_dir, "tresholds.csv"))
        formatter = logging.Formatter(fmt="%(message)s")
        fh.setFormatter(formatter)
        fh.setLevel(logging.INFO)
        logger.handlers = [fh]
        return fh
