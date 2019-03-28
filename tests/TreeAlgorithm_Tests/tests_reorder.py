import unittest
from ddt import ddt
import numpy as np
from tests.context import Pangraph
from tests.context import Node
from tests.context import make_base as n
from tests.context import SequenceID
from tests.context import ConsensusNode
from tests.context import TreePOAConsensusGenerator
from tests.context import MAX1, NODE4
from tests.context import CompatibilityToPath
from tests.context import DataType

@ddt
class TreeAlgorithm_Reconsensus_Test(unittest.TestCase):

    def setUp(self):
        nodes = [Node(node_id=0, base=n('T'), aligned_to=None),
                 Node(node_id=1, base=n('A'), aligned_to=None),
                 Node(node_id=2, base=n('G'), aligned_to=None),
                 Node(node_id=3, base=n('A'), aligned_to=None),
                 Node(node_id=4, base=n('C'), aligned_to=None),
                 Node(node_id=5, base=n('A'), aligned_to=None),
                 Node(node_id=6, base=n('C'), aligned_to=None),

                 Node(node_id=7, base=n('G'), aligned_to=None),
                 Node(node_id=8, base=n('T'), aligned_to=None),
                 Node(node_id=9, base=n('A'), aligned_to=None),
                 Node(node_id=10, base=n('C'), aligned_to=None),
                 Node(node_id=11, base=n('T'), aligned_to=None),
                 Node(node_id=12, base=n('G'), aligned_to=None),
                 Node(node_id=13, base=n('A'), aligned_to=None),
                 Node(node_id=14, base=n('C'), aligned_to=None)
                 ]

        paths = {
            'testseq0': [[0, 1, 2, 3, 4, 5, 6]],
            'testseq1': [[7, 8, 9, 10, 11, 12, 13, 14]]
        }
        self.pangraph = Pangraph(DataType.Nucleotides)
        self.pangraph.nodes = nodes
        self.pangraph.paths = paths

    def test_should_switch_sequences_between_consensuses(self):
        consensus0_path = [0,1,2,3]
        consensus1_path=[7,8,9,10,11]
        original_consensus_nodes = [ConsensusNode(sequences_ids=['testseq1'], consensus_id=0, mincomp=CompatibilityToPath(0,1), consensus_path=consensus0_path),
                              ConsensusNode(sequences_ids=['testseq0'], consensus_id=1, mincomp=CompatibilityToPath(0,1), consensus_path=consensus1_path)]

        tree_generator = TreePOAConsensusGenerator(MAX1([0,1]),NODE4(),"",stop=0.99, re_consensus=True, p=1)

        tree_generator.pangraph = self.pangraph
        expected_reordered = [ConsensusNode(sequences_ids=['testseq0'], consensus_id=0, mincomp=CompatibilityToPath(0.571428571,1), consensus_path=consensus0_path),
                              ConsensusNode(sequences_ids=['testseq1'], consensus_id=1, mincomp=CompatibilityToPath(0.625,1), consensus_path=consensus1_path)]

        actual_reorderd = tree_generator.reorder_consensuses(original_consensus_nodes)
        self.assertEqual(expected_reordered, actual_reorderd)


    def test_should_remove_one_consensus_node(self):
        consensus0_path = [0,1,2,3,4, 5, 7,8,9]
        consensus1_path=[]
        original_consensus_nodes = [ConsensusNode(sequences_ids=['testseq0'], consensus_id=0, mincomp=CompatibilityToPath(0.857142857, 1), consensus_path=consensus0_path),
                                    ConsensusNode(sequences_ids=['testseq1'], consensus_id=1, mincomp=CompatibilityToPath(0,1), consensus_path=consensus1_path)]

        tree_generator = TreePOAConsensusGenerator(MAX1([0,1]),NODE4(),"",stop=0.99, re_consensus=True, p=1)

        tree_generator.pangraph = self.pangraph
        expected_reordered = [ConsensusNode(sequences_ids=['testseq0', 'testseq1'], consensus_id=0, mincomp=CompatibilityToPath(0.375,1), consensus_path=consensus0_path)]

        actual_reorderd = tree_generator.reorder_consensuses(original_consensus_nodes)
        self.assertEqual(expected_reordered, actual_reorderd)

