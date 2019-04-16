import unittest
from ddt import ddt
from tests.context import poa, pNode, pSeq, pPoagraph


def nid(x): return pNode.NodeID(x)


def b(x): return pNode.Base(x)


@ddt
class PangraphPOTranslator_read_top_consensus_Test(unittest.TestCase):

    def setUp(self):
        nodes = [pNode.Node(node_id=nid(0), base=b('T'), aligned_to=None, ),
                 pNode.Node(node_id=nid(1), base=b('A'), aligned_to=nid(2)),
                 pNode.Node(node_id=nid(2), base=b('G'), aligned_to=nid(1)),
                 pNode.Node(node_id=nid(3), base=b('A'), aligned_to=nid(4)),
                 pNode.Node(node_id=nid(4), base=b('C'), aligned_to=nid(3)),
                 pNode.Node(node_id=nid(5), base=b('A'), aligned_to=nid(6)),
                 pNode.Node(node_id=nid(6), base=b('C'), aligned_to=nid(7)),
                 pNode.Node(node_id=nid(7), base=b('G'), aligned_to=nid(8)),
                 pNode.Node(node_id=nid(8), base=b('T'), aligned_to=nid(5)),
                 pNode.Node(node_id=nid(9), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(10), base=b('C'), aligned_to=nid(11)),
                 pNode.Node(node_id=nid(11), base=b('T'), aligned_to=nid(10)),
                 pNode.Node(node_id=nid(12), base=b('G'), aligned_to=None),
                 pNode.Node(node_id=nid(13), base=b('A'), aligned_to=nid(14)),
                 pNode.Node(node_id=nid(14), base=b('C'), aligned_to=nid(13))
                 ]

        sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 3, 5, 9, 10, 13])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [1, 3, 6, 9, 11])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [2, 4, 7, 9, 11, 12])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [2, 4, 8, 9, 11, 12, 14])])],
                              pSeq.SequenceMetadata({'group': '1'})),
        }

        self.poagraph = pPoagraph.Poagraph(nodes, sequences)

    def test_read_consensus_path_seq1_only_in_input(self):
        translator = poa.PoagraphPOTranslator(self.poagraph, [pSeq.SequenceID('seq1')])
        _ = translator.get_input_po_content()

        poa_lines = ["VERSION=pangenome\n",
                     "NAME=pangenome\n",
                     "TITLE=pangenome\n",
                     "LENGTH=5\n",
                     "SOURCECOUNT=2\n",
                     "SOURCENAME=seq1\n",
                     "SOURCEINFO=5 0 100 0 seq1\n",
                     "SOURCENAME=CONSENS0\n",
                     "SOURCEINFO=5 0 100 0 CONSENS0\n",
                     "a:S0S1\n",
                     "a:L0S0S1\n",
                     "c:L1S0S1\n",
                     "a:L2S0S1\n",
                     "t:L2S0S1"]
        actual_consensus_path = translator.read_consensus_paths(poa_lines, [0])
        expected_consensus_path = [1, 3, 6, 9, 11]
        self.assertEqual(expected_consensus_path, actual_consensus_path[0].path)


    def test_subpangraph_construction_from_pangraph_keep_seq_0_1(self):
        translator = poa.PoagraphPOTranslator(self.poagraph, [pSeq.SequenceID('seq0'), pSeq.SequenceID('seq1')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content =  "VERSION=pangenome\n"\
                               "NAME=pangenome\n"\
                               "TITLE=pangenome\n"\
                               "LENGTH=9\n"\
                               "SOURCECOUNT=2\n"\
                               "SOURCENAME=seq0\n"\
                               "SOURCEINFO=7 0 0 -1 seq0\n"\
                               "SOURCENAME=seq1\n"\
                               "SOURCEINFO=5 1 100 -1 seq1\n"\
                               "t:S0\n"\
                               "a:L0S0S1\n"\
                               "a:L1S0S1\n"\
                               "a:L2S0A4\n"\
                               "c:L2S1A3\n"\
                               "a:L3L4S0S1\n"\
                               "c:L5S0A7\n"\
                               "t:L5S1A6\n"\
                               "a:L6S0"
        self.assertEqual(expected_po_content, actual_po_content)

    def test_subpangraph_construction_from_pangraph_keep_seq3(self):
        translator = poa.PoagraphPOTranslator(self.poagraph, [pSeq.SequenceID('seq3')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n" \
                              "NAME=pangenome\n" \
                              "TITLE=pangenome\n" \
                              "LENGTH=7\n" \
                              "SOURCECOUNT=1\n" \
                              "SOURCENAME=seq3\n" \
                              "SOURCEINFO=7 0 100 -1 seq3\n" \
                              "g:S0\n" \
                              "c:L0S0\n" \
                              "t:L1S0\n" \
                              "a:L2S0\n" \
                              "t:L3S0\n" \
                              "g:L4S0\n" \
                              "c:L5S0"
        self.assertEqual(expected_po_content, actual_po_content)

    def test_subpangraph_construction_full_graph(self):
        nodes = [pNode.Node(node_id=nid(0), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(1), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(2), base=b('C'), aligned_to=None),
                 pNode.Node(node_id=nid(3), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(4), base=b('T'), aligned_to=None)
                 ]

        sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4])])],
                              pSeq.SequenceMetadata({'group': '1'}))
        }
        poagraph = pPoagraph.Poagraph(nodes, sequences)
        translator = poa.PoagraphPOTranslator(poagraph, [pSeq.SequenceID('seq0')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n" \
                              "NAME=pangenome\n" \
                              "TITLE=pangenome\n" \
                              "LENGTH=5\n" \
                              "SOURCECOUNT=1\n" \
                              "SOURCENAME=seq0\n" \
                              "SOURCEINFO=5 0 100 -1 seq0\n" \
                              "a:S0\n" \
                              "a:L0S0\n" \
                              "c:L1S0\n" \
                              "a:L2S0\n" \
                              "t:L3S0"
        self.assertEqual(expected_po_content, actual_po_content)

    def test_subpangraph_should_omit_edges_1(self):
        nodes = [pNode.Node(node_id=nid(0), base=b('A'), aligned_to=None),
                          pNode.Node(node_id=nid(1), base=b('C'), aligned_to=None),
                          pNode.Node(node_id=nid(2), base=b('C'), aligned_to=None)]

        sequences = {
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2])])],
                              pSeq.SequenceMetadata({'group': '1'}))
        }
        poagraph = pPoagraph.Poagraph(nodes, sequences)

        translator = poa.PoagraphPOTranslator(poagraph, [pSeq.SequenceID('seq2')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n" \
                              "NAME=pangenome\n" \
                              "TITLE=pangenome\n" \
                              "LENGTH=3\n" \
                              "SOURCECOUNT=1\n" \
                              "SOURCENAME=seq2\n" \
                              "SOURCEINFO=3 0 100 -1 seq2\n" \
                              "a:S0\n" \
                              "c:L0S0\n" \
                              "c:L1S0"
        self.assertEqual(expected_po_content, actual_po_content)


    def test_subpangraph_should_omit_edges_2(self):
        nodes = [pNode.Node(node_id=nid(0), base=b('A'), aligned_to=None),
                          pNode.Node(node_id=nid(1), base=b('C'), aligned_to=None),
                          pNode.Node(node_id=nid(2), base=b('C'), aligned_to=None)]

        sequences = {
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2])])],
                              pSeq.SequenceMetadata({'group': '1'}))
        }
        poagraph = pPoagraph.Poagraph(nodes, sequences)

        translator = poa.PoagraphPOTranslator(poagraph, [pSeq.SequenceID('seq1')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n" \
                              "NAME=pangenome\n" \
                              "TITLE=pangenome\n" \
                              "LENGTH=2\n" \
                              "SOURCECOUNT=1\n" \
                              "SOURCENAME=seq1\n" \
                              "SOURCEINFO=2 0 100 -1 seq1\n" \
                              "a:S0\n" \
                              "c:L0S0"

        self.assertEqual(expected_po_content, actual_po_content)

    def test_subpangraph_should_omit_in_nodes_and_aligned_nodes(self):
        #original pangraph
        nodes = [pNode.Node(node_id=nid(0), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(1), base=b('C'), aligned_to=nid(2)),
                 pNode.Node(node_id=nid(2), base=b('T'), aligned_to=nid(1)),
                 pNode.Node(node_id=nid(3), base=b('G'), aligned_to=None)]

        sequences = {
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 3])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 3])])],
                              pSeq.SequenceMetadata({'group': '1'}))
        }
        poagraph = pPoagraph.Poagraph(nodes, sequences)

        translator = poa.PoagraphPOTranslator(poagraph, [pSeq.SequenceID('seq2')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n" \
                              "NAME=pangenome\n" \
                              "TITLE=pangenome\n" \
                              "LENGTH=3\n" \
                              "SOURCECOUNT=1\n" \
                              "SOURCENAME=seq2\n" \
                              "SOURCEINFO=3 0 100 -1 seq2\n" \
                              "a:S0\n" \
                              "t:L0S0\n" \
                              "g:L1S0"

        self.assertEqual(expected_po_content, actual_po_content)

    def test_subpangraph_unfilled_nodes(self):
        symbol_for_uknown = '?'
        nodes = [pNode.Node(node_id=nid(0), base=b('A'), aligned_to=nid(1)),
                 pNode.Node(node_id=nid(1), base=b('C'), aligned_to=nid(0)),
                 pNode.Node(node_id=nid(2), base=b('G'), aligned_to=None),
                 pNode.Node(node_id=nid(3), base=b(symbol_for_uknown), aligned_to=None),
                 pNode.Node(node_id=nid(4), base=b(symbol_for_uknown), aligned_to=None),
                 pNode.Node(node_id=nid(5), base=b('G'), aligned_to=None),
                 pNode.Node(node_id=nid(6), base=b('C'), aligned_to=None),
                 pNode.Node(node_id=nid(7), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(5), base=b('T'), aligned_to=None)]

        sequences = {
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 3, 4, 7, 8])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [1, 2, 5, 6, 7, 8])])],
                              pSeq.SequenceMetadata({'group': '1'}))
        }
        poagraph = pPoagraph.Poagraph(nodes, sequences)

        translator = poa.PoagraphPOTranslator(poagraph, [pSeq.SequenceID('seq1'),
                                                         pSeq.SequenceID('seq2')])
        actual_po_content = translator.get_input_po_content()
        expected_po_content = "VERSION=pangenome\n" \
                              "NAME=pangenome\n" \
                              "TITLE=pangenome\n" \
                              "LENGTH=9\n" \
                              "SOURCECOUNT=2\n" \
                              "SOURCENAME=seq1\n" \
                              "SOURCEINFO=6 0 100 -1 seq1\n" \
                              "SOURCENAME=seq2\n" \
                              "SOURCEINFO=6 1 100 -1 seq2\n" \
                              "a:S0A1\n" \
                              "c:S1A0\n" \
                              "g:L0L1S0S1\n" \
                              f"{symbol_for_uknown}:L2S0\n" \
                              f"{symbol_for_uknown}:L3S0\n" \
                              "g:L2S1\n" \
                              "c:L5S1\n" \
                              "a:L4L6S0S1\n" \
                              "t:L7S0S1"
        self.assertEqual(expected_po_content, actual_po_content)
