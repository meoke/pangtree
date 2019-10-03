import unittest
from pathlib import Path

from tests.context import pPoagraph, pNode, pSeq, PangenomeFASTA, at_structure, P, pathtools


def nid(x): return pNode.NodeID(x)


def bid(x): return pNode.Base(x)


class ToFASTATests(unittest.TestCase):

    def setUp(self):
        self.fasta_dir = 'tests/output/fasta_files/'

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
                          pNode.Node(node_id=nid(17), base=bid('G'), aligned_to=nid(15))]

        poagraph_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SeqPath([*map(nid, [0, 2, 4, 6, 7, 8, 12, 14, 16])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SeqPath([*map(nid, [3, 4, 6, 7, 10, 12])]),
                               pSeq.SeqPath([*map(nid, [14, 17])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SeqPath([*map(nid, [11])]),
                               pSeq.SeqPath([*map(nid, [13, 14, 15])])],
                              pSeq.SequenceMetadata({'group': '1'})),
        }

        self.poagraph = pPoagraph.Poagraph(poagraph_nodes, poagraph_sequences)

    def test_1_sequences_fasta(self):
        expected_sequences_fasta_path = Path(self.fasta_dir + "_sequences.fasta")

        actual_sequences_fasta_content = PangenomeFASTA.poagraph_to_fasta(self.poagraph)
        expected_sequences_fasta_content = pathtools.get_file_content(expected_sequences_fasta_path)
        self.assertEqual(expected_sequences_fasta_content, actual_sequences_fasta_content)

    def test_2_consensuses_tree_fasta(self):
        expected_consensuses_fasta_path = Path(self.fasta_dir + "consensuses.fasta")

        affinity_tree = at_structure.AffinityTree()
        affinity_tree.nodes = [
            # all members set
            at_structure.AffinityNode(id_=at_structure.AffinityNodeID(0),
                            parent=at_structure.AffinityNodeID(-1),
                            children=[at_structure.AffinityNodeID(1), at_structure.AffinityNodeID(2)],
                            sequences=[pSeq.SequenceID('seq0'),
                                       pSeq.SequenceID('seq1'),
                                       pSeq.SequenceID('seq2'),
                                       pSeq.SequenceID('seq3')],
                            mincomp=at_structure.Compatibility(0.5, P(1)),
                            compatibilities={pSeq.SequenceID('seq0'): at_structure.Compatibility(1.0, P(1)),
                                             pSeq.SequenceID('seq1'): at_structure.Compatibility(0.9, P(1)),
                                             pSeq.SequenceID('seq2'): at_structure.Compatibility(0.95, P(1)),
                                             pSeq.SequenceID('seq3'): at_structure.Compatibility(0.6, P(1))},
                            consensus=pSeq.SeqPath([nid(0), nid(2), nid(5), nid(6),
                                                 nid(10), nid(12), nid(13), nid(16)])),
            # no compatibilities to all, no mincomp
            at_structure.AffinityNode(id_=at_structure.AffinityNodeID(1),
                            parent=at_structure.AffinityNodeID(0),
                            sequences=[pSeq.SequenceID('seq0'),
                                       pSeq.SequenceID('seq1'),
                                       pSeq.SequenceID('seq2')],
                            consensus=pSeq.SeqPath(
                                 [nid(0), nid(2), nid(3), nid(6), nid(10), nid(11), nid(13), nid(17)]))
        ]

        actual_consensuses_fasta_content = PangenomeFASTA.affinity_tree_to_fasta(self.poagraph, affinity_tree)
        expected_consensuses_fasta_content = pathtools.get_file_content(expected_consensuses_fasta_path)
        self.assertEqual(expected_consensuses_fasta_content, actual_consensuses_fasta_content)


if __name__ == '__main__':
    unittest.main()
