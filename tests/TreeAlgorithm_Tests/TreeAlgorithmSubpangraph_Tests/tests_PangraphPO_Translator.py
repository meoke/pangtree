import unittest
from ddt import ddt
import numpy as np
from tests.context import PangraphPO_Translator
from tests.context import Pangraph
from tests.context import Node
from tests.context import nucleotides as n
from tests.context import NodePO, SequencePO

@ddt
class PangraphPO_TranslatorTest(unittest.TestCase):

    def setUp(self):
        nodes = [Node(id=0, base=n.code('T'), aligned_to=None,),
                 Node(id=1, base=n.code('A'), aligned_to=2),
                 Node(id=2, base=n.code('G'), aligned_to=1),
                 Node(id=3, base=n.code('A'), aligned_to=4),
                 Node(id=4, base=n.code('C'), aligned_to=3),
                 Node(id=5, base=n.code('A'), aligned_to=6),
                 Node(id=6, base=n.code('C'), aligned_to=7),
                 Node(id=7, base=n.code('G'), aligned_to=8),
                 Node(id=8, base=n.code('T'), aligned_to=5),
                 Node(id=9, base=n.code('A'), aligned_to=None),
                 Node(id=10, base=n.code('C'), aligned_to=11),
                 Node(id=11, base=n.code('T'), aligned_to=10),
                 Node(id=12, base=n.code('G'), aligned_to=None),
                 Node(id=13, base=n.code('A'), aligned_to=14),
                 Node(id=14, base=n.code('C'), aligned_to=13)
                 ]

        paths = {
            'testseq0': [[0, 1, 3, 5, 9, 10, 13]],
            'testseq1': [[1, 3, 6, 9, 11]],
            'testseq2': [[2, 4, 7, 9, 11, 12]],
            'testseq3': [[2, 4, 8, 9, 11, 12, 14]]
        }
        self.pangraph = Pangraph()
        self.pangraph.nodes = nodes
        self.pangraph.paths = paths

    def test_subpangraph_construction_from_pangraph_keep_seq_0_1(self):
        expected_ponodes = [NodePO(base=4, aligned_to=None, in_nodes=[], sequences_ids=[0]),
                            NodePO(base=1, aligned_to=None, in_nodes=[0], sequences_ids=[0, 1]),
                            NodePO(base=1, in_nodes=[1], aligned_to=None, sequences_ids=[0, 1]),
                            NodePO(base=1, in_nodes=[2], aligned_to=4, sequences_ids=[0]),
                            NodePO(base=2, in_nodes=[2], aligned_to=3, sequences_ids=[1]),
                            NodePO(base=1, in_nodes=[3, 4], aligned_to=None, sequences_ids=[0]),
                            NodePO(base=2, in_nodes=[5], aligned_to=7, sequences_ids=[0]),
                            NodePO(base=4, in_nodes=[5], aligned_to=6, sequences_ids=[1]),
                            NodePO(base=1, in_nodes=[6], aligned_to=None, sequences_ids=[0])]
        expected_posequences = [SequencePO(name='testseq0',
                                           nodes_count=7,
                                           weight=-1,
                                           consensus_id=-1,
                                           start_node_id=0),
                                SequencePO(name='testseq1',
                                           nodes_count=5,
                                           weight=-1,
                                           consensus_id=-1,
                                           start_node_id=1),
                                ]

        expected_nodes_ids_mapping = {0: 0, 1: 1, 2: 3, 3: 5, 4: 6, 5: 9, 6: 10, 7: 11, 8: 13}

        translator = PangraphPO_Translator(self.pangraph, ['testseq0', 'testseq1'])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n"\
                               "NAME=pangenome\n"\
                               "TITLE=pangenome\n"\
                               "LENGTH=9\n"\
                               "SOURCECOUNT=2\n"\
                               "SOURCENAME=testseq0\n"\
                               "SOURCEINFO=7 0 100 -1 testseq0\n"\
                               "SOURCENAME=testseq1\n"\
                               "SOURCEINFO=5 1 0 -1 testseq1\n"\
                               "t:S0\n"\
                               "a:L0S0S1\n"\
                               "a:L1S0S1\n"\
                               "a:L2S0A4\n"\
                               "c:L2S1A3\n"\
                               "a:L3L4S0\n"\
                               "c:L5S0A7\n"\
                               "t:L5S1A6\n"\
                               "a:L6S0\n"
        self.assertEqual(expected_po_content, actual_po_content)

    def test_subpangraph_construction_from_pangraph_keep_seq3(self):
        expected_nodes = [Node(id=0, base=n.code('G'), in_nodes=[], aligned_to=None),
                          Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
                          Node(id=2, base=n.code('T'), in_nodes=[1], aligned_to=None),
                          Node(id=3, base=n.code('A'), in_nodes=[2], aligned_to=None),
                          Node(id=4, base=n.code('T'), in_nodes=[3], aligned_to=None),
                          Node(id=5, base=n.code('G'), in_nodes=[4], aligned_to=None),
                          Node(id=6, base=n.code('C'), in_nodes=[5], aligned_to=None)
                          ]

        expected_paths_to_node_ids = {
            'testseq3': [0, 1, 2, 3, 4, 5, 6]
        }
        expected_pangraph = Pangraph()
        expected_pangraph.update_nodes(expected_nodes)
        expected_pangraph.set_paths(len(expected_nodes), expected_paths_to_node_ids)

        expected_nodes_ids_mapping = {0: 2, 1: 4, 2: 8, 3: 9, 4: 11, 5: 12, 6: 14}

        expected_subpangraph = SubPangraph(Pangraph(), [])
        expected_subpangraph.pangraph = expected_pangraph
        expected_subpangraph.nodes_ids_mapping = expected_nodes_ids_mapping

        actual_subpangraph = SubPangraph(self.pangraph, ['testseq3'])
        self.compare_subpangraphs(expected_subpangraph, actual_subpangraph)


    def test_subsubpangraph_construction(self):
        expected_nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
                          Node(id=1, base=n.code('A'), in_nodes=[0], aligned_to=None),
                          Node(id=2, base=n.code('C'), in_nodes=[1], aligned_to=None),
                          Node(id=3, base=n.code('A'), in_nodes=[2], aligned_to=None),
                          Node(id=4, base=n.code('T'), in_nodes=[3], aligned_to=None)
                          ]

        expected_paths_to_node_ids = {
            'testseq1': [0, 1, 2, 3, 4]
        }
        expected_pangraph = Pangraph()
        expected_pangraph.update_nodes(expected_nodes)
        expected_pangraph.set_paths(len(expected_nodes), expected_paths_to_node_ids)

        expected_nodes_ids_mapping = {0: 1, 1: 2, 2: 4, 3: 5, 4: 7}

        expected_subpangraph = SubPangraph(Pangraph(), [])
        expected_subpangraph.pangraph = expected_pangraph
        expected_subpangraph.nodes_ids_mapping = expected_nodes_ids_mapping

        actual_subpangraph = SubPangraph(self.pangraph, ['testseq0', 'testseq1'])
        actual_subpangraph = actual_subpangraph.keep_sources_ids(['testseq1'])  # tu pewnie musi być wg names
        self.compare_subpangraphs(expected_subpangraph, actual_subpangraph)

    def test_subsubpangraph_construction_2(self):
        expected_nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
                          Node(id=1, base=n.code('A'), in_nodes=[0], aligned_to=None),
                          Node(id=2, base=n.code('C'), in_nodes=[1], aligned_to=None),
                          Node(id=3, base=n.code('A'), in_nodes=[2], aligned_to=None),
                          Node(id=4, base=n.code('T'), in_nodes=[3], aligned_to=None)
                          ]

        expected_paths_to_node_ids = {
            'testseq1': [0, 1, 2, 3, 4]
        }
        expected_pangraph = Pangraph()
        expected_pangraph.update_nodes(expected_nodes)
        expected_pangraph.set_paths(len(expected_nodes), expected_paths_to_node_ids)

        expected_nodes_ids_mapping = {0: 1, 1: 2, 2: 4, 3: 5, 4: 7}

        expected_subpangraph = SubPangraph(Pangraph(), [])
        expected_subpangraph.pangraph = expected_pangraph
        expected_subpangraph.nodes_ids_mapping = expected_nodes_ids_mapping

        actual_subpangraph = SubPangraph(self.pangraph, ['testseq0', 'testseq1'])
        actual_subpangraph = actual_subpangraph.keep_sources_ids(['testseq1'])  # tu pewnie musi być wg names
        self.compare_subpangraphs(expected_subpangraph, actual_subpangraph)

    def test_subpangraph_should_omit_edges_1(self):
        pangraph_nodes = [Node(id=0,base=n.code('A'), in_nodes=[], aligned_to=None),
                          Node(id=1,base=n.code('C'), in_nodes=[0], aligned_to=None),
                          Node(id=2,base=n.code('C'), in_nodes=[0, 1], aligned_to=None)]
        pangraph_paths_to_nodes_ids = {
            'seq1': [0, 2],
            'seq2': [0, 1, 2]
        }
        pangraph = Pangraph()
        pangraph.update_nodes(pangraph_nodes)
        pangraph.set_paths(len(pangraph_nodes), pangraph_paths_to_nodes_ids)

        # remove seq1
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
            Node(id=2, base=n.code('C'), in_nodes=[1], aligned_to=None)
        ]
        expected_paths_to_nodes_ids = {
            'seq2': [0, 1, 2]
        }
        expected_nodes_ids_mapping = {0: 0, 1: 1, 2: 2}

        expected_pangraph = Pangraph()
        expected_pangraph.update_nodes(expected_nodes)
        expected_pangraph.set_paths(len(expected_nodes), expected_paths_to_nodes_ids)

        expected_subpangraph = SubPangraph(Pangraph(), [])
        expected_subpangraph.pangraph = expected_pangraph
        expected_subpangraph.nodes_ids_mapping = expected_nodes_ids_mapping

        actual_subpangraph = SubPangraph(pangraph, ['seq2'])
        self.compare_subpangraphs(expected_subpangraph, actual_subpangraph)


    def test_subpangraph_should_omit_edges_2(self):
        pangraph_nodes = [Node(id=0,base=n.code('A'), in_nodes=[], aligned_to=None),
                          Node(id=1,base=n.code('C'), in_nodes=[0], aligned_to=None),
                          Node(id=2,base=n.code('C'), in_nodes=[0, 1], aligned_to=None)]
        pangraph_paths_to_nodes_ids = {
            'seq1': [0, 2],
            'seq2': [0, 1, 2]
        }
        pangraph = Pangraph()
        pangraph.update_nodes(pangraph_nodes)
        pangraph.set_paths(len(pangraph_nodes), pangraph_paths_to_nodes_ids)

        # remove seq2
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
        ]
        expected_paths_to_nodes_ids = {
            'seq1': [0, 1]
        }
        expected_nodes_ids_mapping = {0: 0, 1: 2}

        expected_pangraph = Pangraph()
        expected_pangraph.update_nodes(expected_nodes)
        expected_pangraph.set_paths(len(expected_nodes), expected_paths_to_nodes_ids)

        expected_subpangraph = SubPangraph(Pangraph(), [])
        expected_subpangraph.pangraph = expected_pangraph
        expected_subpangraph.nodes_ids_mapping = expected_nodes_ids_mapping

        actual_subpangraph = SubPangraph(pangraph, ['seq1'])
        self.compare_subpangraphs(expected_subpangraph, actual_subpangraph)

    def test_subpangraph_should_omit_in_nodes_and_aligned_nodes(self):
        #original pangraph
        pangraph_nodes = [Node(id=0,base=n.code('A'), in_nodes=[], aligned_to=None),
                          Node(id=1,base=n.code('C'), in_nodes=[0], aligned_to=2),
                          Node(id=2,base=n.code('T'), in_nodes=[0], aligned_to=1),
                          Node(id=3,base=n.code('G'), in_nodes=[1, 2], aligned_to=None)]
        pangraph_paths_to_nodes_ids = {
            'seq1': [0, 1, 3],
            'seq2': [0, 2, 3]
        }
        pangraph = Pangraph()
        pangraph.update_nodes(pangraph_nodes)
        pangraph.set_paths(len(pangraph_nodes), pangraph_paths_to_nodes_ids)

        # remove seq1
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('T'), in_nodes=[0], aligned_to=None),
            Node(id=2, base=n.code('G'), in_nodes=[1], aligned_to=None)
        ]
        expected_paths_to_nodes_ids = {
            'seq2': [0, 1, 2]
        }
        expected_nodes_ids_mapping = {0: 0, 1: 2, 2: 3}

        expected_pangraph = Pangraph()
        expected_pangraph.update_nodes(expected_nodes)
        expected_pangraph.set_paths(len(expected_nodes), expected_paths_to_nodes_ids)

        expected_subpangraph = SubPangraph(Pangraph(), [])
        expected_subpangraph.pangraph = expected_pangraph
        expected_subpangraph.nodes_ids_mapping = expected_nodes_ids_mapping

        actual_subpangraph = SubPangraph(pangraph, ['seq2'])
        self.compare_subpangraphs(expected_subpangraph, actual_subpangraph)


    def compare_subpangraphs(self, expected_subpangraph, actual_subpangraph):
        self.assertEqual(expected_subpangraph.pangraph._nodes, actual_subpangraph.pangraph._nodes)
        self.assertTrue(np.array_equal(expected_subpangraph.pangraph._pathmanager.paths, actual_subpangraph.pangraph._pathmanager.paths))
        self.assertEqual(expected_subpangraph.nodes_ids_mapping, actual_subpangraph.nodes_ids_mapping)