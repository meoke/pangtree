import unittest
from pathlib import Path


from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import missings, msa
from pangtreebuild.tools import pathtools
from pangtreebuild.serialization import fasta
from pangtreebuild.affinity_tree import tree
from pangtreebuild.affinity_tree import parameters as at_params

def nid(x): return graph.NodeID(x)


def bid(x): return graph.Base(x)


class ToFASTATests(unittest.TestCase):

    def setUp(self):
        self.fasta_dir = Path(__file__).parent.joinpath('fasta_files/')

        poagraph_nodes = [graph.Node(node_id=nid(0), base=bid('A'), aligned_to=nid(1)),
                          graph.Node(node_id=nid(1), base=bid('G'), aligned_to=nid(0)),
                          graph.Node(node_id=nid(2), base=bid('C'), aligned_to=nid(3)),
                          graph.Node(node_id=nid(3), base=bid('G'), aligned_to=nid(2)),
                          graph.Node(node_id=nid(4), base=bid('A'), aligned_to=nid(5)),
                          graph.Node(node_id=nid(5), base=bid('T'), aligned_to=nid(4)),
                          graph.Node(node_id=nid(6), base=bid('G'), aligned_to=None),
                          graph.Node(node_id=nid(7), base=bid('G'), aligned_to=None),
                          graph.Node(node_id=nid(8), base=bid('A'), aligned_to=nid(9)),
                          graph.Node(node_id=nid(9), base=bid('C'), aligned_to=nid(10)),
                          graph.Node(node_id=nid(10), base=bid('G'), aligned_to=nid(11)),
                          graph.Node(node_id=nid(11), base=bid('T'), aligned_to=nid(8)),
                          graph.Node(node_id=nid(12), base=bid('A'), aligned_to=nid(13)),
                          graph.Node(node_id=nid(13), base=bid('C'), aligned_to=nid(12)),
                          graph.Node(node_id=nid(14), base=bid('T'), aligned_to=None),
                          graph.Node(node_id=nid(15), base=bid('A'), aligned_to=nid(16)),
                          graph.Node(node_id=nid(16), base=bid('C'), aligned_to=nid(17)),
                          graph.Node(node_id=nid(17), base=bid('G'), aligned_to=nid(15))]

        poagraph_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0, 2, 4, 6, 7, 8, 12, 14, 16])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [3, 4, 6, 7, 10, 12])]),
                                graph.SeqPath([*map(nid, [14, 17])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [11])]),
                                graph.SeqPath([*map(nid, [13, 14, 15])])],
                               graph.SequenceMetadata({'group': '1'})),
        }

        self.poagraph = graph.Poagraph(poagraph_nodes, poagraph_sequences)

    def test_1_sequences_fasta(self):
        expected_sequences_fasta_path = self.fasta_dir.joinpath("sequences.fasta")

        actual_sequences_fasta_content = fasta.poagraph_to_fasta(self.poagraph)
        expected_sequences_fasta_content = pathtools.get_file_content(expected_sequences_fasta_path)
        self.assertEqual(expected_sequences_fasta_content, actual_sequences_fasta_content)

    def test_2_consensuses_tree_fasta(self):
        expected_consensuses_fasta_path = self.fasta_dir.joinpath("consensuses.fasta")

        affinity_tree = tree.AffinityTree()
        affinity_tree.nodes = [
            # all members set
            tree.AffinityNode(id_=tree.AffinityNodeID(0),
                              parent=tree.AffinityNodeID(-1),
                              children=[tree.AffinityNodeID(1), tree.AffinityNodeID(2)],
                              sequences=[msa.SequenceID('seq0'),
                                         msa.SequenceID('seq1'),
                                         msa.SequenceID('seq2'),
                                         msa.SequenceID('seq3')],
                              mincomp=graph.Compatibility(0.5, at_params.P(1)),
                              compatibilities={msa.SequenceID('seq0'): graph.Compatibility(1.0, at_params.P(1)),
                                               msa.SequenceID('seq1'): graph.Compatibility(0.9, at_params.P(1)),
                                               msa.SequenceID('seq2'): graph.Compatibility(0.95, at_params.P(1)),
                                               msa.SequenceID('seq3'): graph.Compatibility(0.6, at_params.P(1))},
                              consensus=graph.SeqPath([nid(0), nid(2), nid(5), nid(6),
                                                       nid(10), nid(12), nid(13), nid(16)])),
            # no compatibilities to all, no mincomp
            tree.AffinityNode(id_=tree.AffinityNodeID(1),
                              parent=tree.AffinityNodeID(0),
                              sequences=[msa.SequenceID('seq0'),
                                         msa.SequenceID('seq1'),
                                         msa.SequenceID('seq2')],
                              consensus=graph.SeqPath([nid(0), nid(2), nid(3),
                                                       nid(6), nid(10), nid(11),
                                                       nid(13), nid(17)]))
        ]

        actual_consensuses_fasta_content = fasta.affinity_tree_to_fasta(self.poagraph, affinity_tree)
        expected_consensuses_fasta_content = pathtools.get_file_content(expected_consensuses_fasta_path)
        self.assertEqual(expected_consensuses_fasta_content, actual_consensuses_fasta_content)


if __name__ == '__main__':
    unittest.main()
