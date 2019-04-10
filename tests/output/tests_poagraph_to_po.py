import unittest
from pathlib import Path

from tests.context import pathtools, pNode, pSeq, pPoagraph, PangenomePO


def nid(x): return pNode.NodeID(x)


def bid(x): return pNode.Base(x)


class PoagraphToPOTests(unittest.TestCase):

    def setUp(self):
        self.po_files_dir = 'tests/output/po_files/'

    def test_1_typical_pangraph(self):
        expected_po_content_path = Path(self.po_files_dir + "test_1.po")

        poagraph_nodes = [pNode.Node(node_id=nid(0), base=bid('A'), aligned_to=nid(1)),
                          pNode.Node(node_id=nid(1), base=bid('G'), aligned_to=nid(0)),
                          pNode.Node(node_id=nid(2), base=bid('C'), aligned_to=nid(3)),
                          pNode.Node(node_id=nid(3), base=bid('G'), aligned_to=nid(2)),
                          pNode.Node(node_id=nid(4), base=bid('A'), aligned_to=nid(5)),
                          pNode.Node(node_id=nid(5), base=bid('T'), aligned_to=nid(4)),
                          pNode.Node(node_id=nid(6), base=bid('G'), aligned_to=None),
                          pNode.Node(node_id=nid(7), base=bid('G'), aligned_to=None),
                          pNode.Node(node_id=nid(8), base=bid('A'), aligned_to=nid(9)),
                          pNode.Node(node_id=nid(9), base=bid('C'), aligned_to=nid(10)),
                          pNode.Node(node_id=nid(10), base=bid('G'), aligned_to=nid(11)),
                          pNode.Node(node_id=nid(11), base=bid('T'), aligned_to=nid(8)),
                          pNode.Node(node_id=nid(12), base=bid('A'), aligned_to=nid(13)),
                          pNode.Node(node_id=nid(13), base=bid('C'), aligned_to=nid(12)),
                          pNode.Node(node_id=nid(14), base=bid('T'), aligned_to=None),
                          pNode.Node(node_id=nid(15), base=bid('A'), aligned_to=nid(16)),
                          pNode.Node(node_id=nid(16), base=bid('C'), aligned_to=nid(17)),
                          pNode.Node(node_id=nid(17), base=bid('G'), aligned_to=nid(15))
                         ]

        poagraph_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 4, 6, 7, 8, 12, 14, 16])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [1, 2, 5, 6, 7, 9])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [3, 4, 6, 7, 10, 12, 14, 17])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [11, 13, 14, 15])])],
                              pSeq.SequenceMetadata({'group': '1'})),
        }

        poagraph = pPoagraph.Poagraph(poagraph_nodes, poagraph_sequences)

        actual_po_content = PangenomePO.poagraph_to_PangenomePO(poagraph)
        expected_po_content = pathtools.get_file_content(expected_po_content_path)
        self.assertEqual(expected_po_content, actual_po_content)

    def test_2_consensuses_and_empty_sequences(self):
        expected_po_content_path = Path(self.po_files_dir + "test_2.po")

        poagraph_nodes = [pNode.Node(node_id=nid(0), base=bid('C'), aligned_to=nid(1)),
                          pNode.Node(node_id=nid(1), base=bid('T'), aligned_to=nid(0)),
                          pNode.Node(node_id=nid(2), base=bid('A'), aligned_to=nid(3)),
                          pNode.Node(node_id=nid(3), base=bid('G'), aligned_to=nid(2)),
                          pNode.Node(node_id=nid(4), base=bid('C'), aligned_to=None),
                          pNode.Node(node_id=nid(5), base=bid('T'), aligned_to=None),
                          pNode.Node(node_id=nid(6), base=bid('A'), aligned_to=nid(7)),
                          pNode.Node(node_id=nid(7), base=bid('T'), aligned_to=nid(6)),
                          pNode.Node(node_id=nid(8), base=bid('G'), aligned_to=None)
                          ]

        poagraph_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 3, 4, 5, 6, 8])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [1, 2, 4, 5, 7, 8])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('CONSENS0'):
                pSeq.Sequence(pSeq.SequenceID('CONSENS0'),
                              [pSeq.SequencePath([*map(nid, [0, 3, 4, 5, 7, 8])])],
                              None),
            pSeq.SequenceID('CONSENS1'):
                pSeq.Sequence(pSeq.SequenceID('CONSENS1'),
                              [pSeq.SequencePath([*map(nid, [1, 2, 4, 5, 6, 8])])],
                              None),
        }

        poagraph = pPoagraph.Poagraph(poagraph_nodes, poagraph_sequences)

        actual_po_content = PangenomePO.poagraph_to_PangenomePO(poagraph)
        expected_po_content = pathtools.get_file_content(expected_po_content_path)
        self.assertEqual(expected_po_content, actual_po_content)



if __name__ == '__main__':
    unittest.main()
