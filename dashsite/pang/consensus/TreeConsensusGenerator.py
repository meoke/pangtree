from pathlib import Path
from collections import deque
from typing import List

from consensus.ConsensusesTree import ConsensusesTree
from consensus.ConsensusNode import ConsensusNode
from consensus.algorithm.FindCutoff import FindMaxCutoff, FindNodeCutoff
from graph2.Pangraph import Pangraph
from graph2.custom_types import SequenceID, NodeID
from metadata.MultialignmentMetadata import MultialignmentMetadata
import consensus.poapoa as poapoa
from consensus.errors import TreeConsensusGenerationException


class TreeConsensusGenerator:
    def __init__(self,
                 max_node_strategy: FindMaxCutoff,
                 node_cutoff_strategy: FindNodeCutoff,
                 stop: float,
                 re_consensus):
        self.max_node_strategy = max_node_strategy
        self.node_cutoff_strategy = node_cutoff_strategy
        self.stop = stop
        self.re_consensus = re_consensus
        self.pangraph = None
        self.genomes_info = None
        self.output_dir = None
        self.consensuses_tree = None

    def get_consensuses_tree(self,
                             pangraph: Pangraph,
                             genomes_info: MultialignmentMetadata,
                             output_dir: Path
                             ) -> ConsensusesTree:
        self.pangraph = pangraph
        self.genomes_info = genomes_info
        self.output_dir = output_dir

        if not self.pangraph_is_valid():
            raise TreeConsensusGenerationException("Invalid pangraph. Cannot find consensuses.")

        self.consensuses_tree = self.init_consensuses_tree()

        nodes_to_process = deque([self.consensuses_tree.get_node(0)])
        while nodes_to_process:
            node = nodes_to_process.pop()
            children_nodes = self.get_children_nodes(node)

            if len(children_nodes) == 1:
                continue

            for child in children_nodes:
                child.compatibilities_to_all = pangraph.get_compatibilities(sequences_ids=list(pangraph.paths.keys()),
                                                                            consensus=child.consensus_path)
                node.children_nodes_ids.append(child.consensus_id)
                self.consensuses_tree.nodes.append(child)
                if not self.node_is_ready(child):
                    nodes_to_process.append(child)
        return self.consensuses_tree

    def init_consensuses_tree(self):
        consensuses_tree = ConsensusesTree()
        root_node = self.get_root_node()
        consensuses_tree.nodes.append(root_node)
        return consensuses_tree

    def get_top_consensus(self, sequences_ids: List[SequenceID], name: str):
        try:
            return poapoa.get_top_consensus(pangraph=self.pangraph,
                                                      sequences_ids=sequences_ids,
                                                      output_dir=self.output_dir,
                                                      file_prefix=name)
        except TreeConsensusGenerationException as e:
            raise TreeConsensusGenerationException(f'Cannot find {name} consensus.') from e

    def get_root_node(self) -> ConsensusNode:
        consensus_path = self.get_top_consensus(sequences_ids=self.pangraph.paths.keys(), name="root")
        compatibilities = self.pangraph.get_compatibilities(sequences_ids=self.pangraph.paths.keys(),
                                                            consensus=consensus_path)
        return ConsensusNode(parent_node_id=None,
                             sequences_ids=[*self.pangraph.paths.keys()],
                             consensus_id=0,
                             mincomp=self.get_mincomp(self.pangraph.paths.keys(), compatibilities),
                             compatibilities_to_all=compatibilities,
                             consensus_path=consensus_path)

    def pangraph_is_valid(self):
        if len(self.pangraph.paths) > 0:
            return True
        return False

    def node_is_ready(self, node: ConsensusNode) -> bool:
        if len(node.sequences_ids) == 1 or node.mincomp >= self.stop:
            return True
        return False

    def get_children_nodes(self, node: ConsensusNode):
        children_nodes: List[ConsensusNode] = []

        current_path_names = set(node.sequences_ids)

        so_far_cutoffs = []
        while current_path_names:
            consensus_path = self.get_top_consensus(sequences_ids=current_path_names,
                                                    name=f"{node.consensus_id}_{len(current_path_names)}_all")
            compatibilities_to_consensus = self.pangraph.get_compatibilities(sequences_ids=current_path_names,
                                                                             consensus=consensus_path)

            max_sequences_ids = self.get_max_compatible_sequences_ids(compatibilities_to_consensus)
            max_consensus_path = self.get_top_consensus(sequences_ids=max_sequences_ids,
                                                        name=f"{node.consensus_id}_{len(current_path_names)}_max")
            compatibilities_to_max_consensus = self.pangraph.get_compatibilities(sequences_ids=current_path_names,
                                                                                 consensus=max_consensus_path)

            qualified_sequences_ids, node_cutoff = self.get_qualified_sequences_ids_and_cutoff(
                                                        compatibilities_to_max_consensus=compatibilities_to_max_consensus,
                                                        so_far_cutoffs=so_far_cutoffs)


            consensus_node = ConsensusNode(parent_node_id=node.consensus_id,
                                           sequences_ids=qualified_sequences_ids,
                                           consensus_id=len(self.consensuses_tree.nodes) + len(so_far_cutoffs),
                                           mincomp=self.get_mincomp(compatibtle_sequences_ids=qualified_sequences_ids,
                                                                    compatibilities_to_consensus=compatibilities_to_max_consensus),
                                           consensus_path=max_consensus_path)
            so_far_cutoffs.append(node_cutoff)
            children_nodes.append(consensus_node)
            current_path_names = list(set(current_path_names) - set(qualified_sequences_ids))

        if self.re_consensus:
            children_nodes = self.reorder_consensuses(children_nodes)

        return children_nodes

    def get_max_compatible_sequences_ids(self, compatibilities_to_consensus):
        max_cutoff = self.max_node_strategy.find_max_cutoff([*compatibilities_to_consensus.values()], log=True)
        max_sequences_ids = self.get_sequences_ids_above_cutoff(compatibilities_to_consensus, max_cutoff)
        return max_sequences_ids

    def get_qualified_sequences_ids_and_cutoff(self, compatibilities_to_max_consensus, so_far_cutoffs):
        node_cutoff = self.node_cutoff_strategy.find_node_cutoff(compatibilities=[*compatibilities_to_max_consensus.values()],
                                                                 so_far_cutoffs=so_far_cutoffs,
                                                                 log=True)
        compatibtle_sequences_ids = self.get_sequences_ids_above_cutoff(compatibilities_to_max_consensus, node_cutoff)
        return compatibtle_sequences_ids, node_cutoff

    def get_mincomp(self, compatibtle_sequences_ids, compatibilities_to_consensus):
        return min([comp
                for seq_id, comp
                in compatibilities_to_consensus.items()
                if seq_id in compatibtle_sequences_ids])

    def get_sequences_ids_above_cutoff(self, compatibilities, cutoff):
        return [seq_id for seq_id, comp in compatibilities.items() if comp >= cutoff]

    def reorder_consensuses(self, children_nodes):
        return children_nodes
