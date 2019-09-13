import unittest
from pathlib import Path

from tests.context import pPoagraph, pNode, pSeq, PangenomeFASTA, CT, P, pathtools


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
                              [pSeq.SequencePath([*map(nid, [0, 2, 4, 6, 7, 8, 12, 14, 16])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [3, 4, 6, 7, 10, 12])]),
                               pSeq.SequencePath([*map(nid, [14, 17])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [11])]),
                               pSeq.SequencePath([*map(nid, [13, 14, 15])])],
                              pSeq.SequenceMetadata({'group': '1'})),
        }

        self.poagraph = pPoagraph.Poagraph(poagraph_nodes, poagraph_sequences)

    def test_1_sequences_fasta(self):
        expected_sequences_fasta_path = Path(self.fasta_dir + "sequences.fasta")

        actual_sequences_fasta_content = PangenomeFASTA.poagraph_to_fasta(self.poagraph)
        expected_sequences_fasta_content = pathtools.get_file_content(expected_sequences_fasta_path)
        self.assertEqual(expected_sequences_fasta_content, actual_sequences_fasta_content)

    def test_2_consensuses_tree_fasta(self):
        expected_consensuses_fasta_path = Path(self.fasta_dir + "consensuses.fasta")

        consensuses_tree = CT.ConsensusTree()
        consensuses_tree.nodes = [
            # all members set
            CT.ConsensusNode(consensus_id=CT.ConsensusNodeID(0),
                             parent_node_id=CT.ConsensusNodeID(-1),
                             children_nodes_ids=[CT.ConsensusNodeID(1), CT.ConsensusNodeID(2)],
                             sequences_ids=[pSeq.SequenceID('seq0'),
                                            pSeq.SequenceID('seq1'),
                                            pSeq.SequenceID('seq2'),
                                            pSeq.SequenceID('seq3')],
                             mincomp=CT.CompatibilityToPath(0.5, P(1)),
                             compatibilities_to_all={pSeq.SequenceID('seq0'): CT.CompatibilityToPath(1.0, P(1)),
                                                     pSeq.SequenceID('seq1'): CT.CompatibilityToPath(0.9, P(1)),
                                                     pSeq.SequenceID('seq2'): CT.CompatibilityToPath(0.95, P(1)),
                                                     pSeq.SequenceID('seq3'): CT.CompatibilityToPath(0.6, P(1))},
                             consensus_path=pSeq.SequencePath([nid(0), nid(2), nid(5), nid(6),
                                                               nid(10), nid(12), nid(13), nid(16)])),
            # no compatibilities to all, no mincomp
            CT.ConsensusNode(consensus_id=CT.ConsensusNodeID(1),
                             parent_node_id=CT.ConsensusNodeID(0),
                             sequences_ids=[pSeq.SequenceID('seq0'),
                                            pSeq.SequenceID('seq1'),
                                            pSeq.SequenceID('seq2')],
                             consensus_path=pSeq.SequencePath(
                                 [nid(0), nid(2), nid(3), nid(6), nid(10), nid(11), nid(13), nid(17)]))
        ]

        actual_consensuses_fasta_content = PangenomeFASTA.consensuses_tree_to_fasta(self.poagraph, consensuses_tree)
        expected_consensuses_fasta_content = pathtools.get_file_content(expected_consensuses_fasta_path)
        self.assertEqual(expected_consensuses_fasta_content, actual_consensuses_fasta_content)


if __name__ == '__main__':
    unittest.main()
