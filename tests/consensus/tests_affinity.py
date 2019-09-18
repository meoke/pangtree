import unittest
from typing import List

from ddt import unpack, data, ddt

from tests.context import pNode, pSeq, pPoagraph, CT, P, cutoffs, tree_generator


def sid(x): return pSeq.SequenceID(x)


def nid(x): return pNode.NodeID(x)


def b(x): return pNode.Base(x)


@ddt
class AffinityTreeGenerationTests(unittest.TestCase):

    @data((P(0.5), CT.CompatibilityToPath(0.836660026534076)),
          (P(1), CT.CompatibilityToPath(0.7)),
          (P(4), CT.CompatibilityToPath(0.6561)))
    @unpack
    def test_1_p_parameter_influence(self, p: P, expected_cutoff: CT.CompatibilityToPath):
        nodes = [pNode.Node(node_id=nid(0), base=b('T'), aligned_to=None),
                 pNode.Node(node_id=nid(1), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(2), base=b('G'), aligned_to=None),
                 pNode.Node(node_id=nid(3), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(4), base=b('C'), aligned_to=None),
                 pNode.Node(node_id=nid(5), base=b('A'), aligned_to=None),
                 pNode.Node(node_id=nid(6), base=b('C'), aligned_to=None),
                 pNode.Node(node_id=nid(7), base=b('G'), aligned_to=None),
                 pNode.Node(node_id=nid(8), base=b('T'), aligned_to=None),
                 pNode.Node(node_id=nid(9), base=b('A'), aligned_to=None)
                 ]

        sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 18, 9])])],
                              pSeq.SequenceMetadata({})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 8, 9])])],
                              pSeq.SequenceMetadata({})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 7, 8, 9])])],
                              pSeq.SequenceMetadata({})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 3, 4, 5, 6, 7, 8, 9])])],
                              pSeq.SequenceMetadata({})),
            pSeq.SequenceID('seq4'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 2, 3, 4, 5, 6, 7, 8, 9])])],
                              pSeq.SequenceMetadata({}))
        }

        poagraph = pPoagraph.Poagraph(nodes, sequences)

        consensus_path = pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 18, 19])])
        compatibilities = poagraph.get_compatibilities(poagraph.get_sequences_ids(), consensus_path, p)

        actual_cutoff = tree_generator.find_node_cutoff([c for c in compatibilities.values()], []).cutoff
        self.assertAlmostEqual(expected_cutoff.value, actual_cutoff.value)

    @data(
        # single compatibility value
        (0.5, [CT.CompatibilityToPath(0.5)]),

        # two compatibilities values
        (0.7, [CT.CompatibilityToPath(0.5), CT.CompatibilityToPath(0.7)]),
        (1, [CT.CompatibilityToPath(1), CT.CompatibilityToPath(0.45)]),
        (0.9, [CT.CompatibilityToPath(0.9), CT.CompatibilityToPath(0.5)]),

        # repeated values
        (0.7, [*map(CT.CompatibilityToPath, [0.5, 0.7, 0.7])]),
        (0.9, [*map(CT.CompatibilityToPath, [0.9, 0.5, 0.5])]),
        (1, [*map(CT.CompatibilityToPath, [0.45, 1, 0.45, 0.45])]),

        # many unique compatibilities values
        (.8, [*map(CT.CompatibilityToPath, [.3, .4, .8])]),
        (0.91, [*map(CT.CompatibilityToPath, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97])]),
        (0.91, [*map(CT.CompatibilityToPath, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97])]),
        (1, [*map(CT.CompatibilityToPath, [0.81, 0.75, 0.8, 0.81, 1])]),
        (0.9, [*map(CT.CompatibilityToPath, [0.5, 0.9, 0.99])]),
        (0.7, [*map(CT.CompatibilityToPath, [0.2, 0.85, 0.7, 0.8])]),
        (0.99, [*map(CT.CompatibilityToPath, [0.99, 0.9, 0.99])]),
        (0.99, [*map(CT.CompatibilityToPath, [0.99])]),

        # repeated distance between values
        (.4, [*map(CT.CompatibilityToPath, [.3, .4, .5])]),

        # all the same values
        (.1, [*map(CT.CompatibilityToPath, [.1, .1, .1])])
    )
    @unpack
    def test_2_find_cutoff_no_so_far_values(self, expected_cutoff: float, compatibilities: List[CT.CompatibilityToPath]):
        actual_cutoff = tree_generator.find_node_cutoff(compatibilities, []).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)

    def test_3_find_cutoff_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            _ = tree_generator.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")


    @data(
        # guard <= all compatibilities
        (0.2, [0.2, 0.7, 0.8, 0.85], [0.1, 0.01, 0]),
        (0.7, [0.7, 0.85, 0.7, 0.8], [0.1, 0.01, 0]),
        (0.8, [0.7, 0.7, 0.85, 0.8], [0.85, 0.91, 1.0]),

        # guard > all compatibilities
        (0.6, [0.3, 0.6, 0.61, 0.61], [0.99]),  # gap to guard bigger than winning gap
        (0.9, [0.2, 0.97, 0.98, 0.9], [0.99]),  # gap to guard smaller than winning gap

        # guard between compatibilities
        (0.5, [0.2, 0.57, 0.58, 0.5], [0.55]),  # take smaller than guard
        (0.58, [0.2, 0.27, 0.58, 0.2], [0.55]),  # take greater than guard
        (0.55, [0.2, 0.58, 0.27, 0.55], [0.55])  # take equal to guard
    )
    @unpack
    def test_4_find_cutoff_with_so_far_values(self, expected_cutoff, compatibilities, so_far_cutoffs):
        compatibilities = [CT.CompatibilityToPath(c) for c in compatibilities]
        so_far_cutoffs = [CT.CompatibilityToPath(c) for c in so_far_cutoffs]
        actual_cutoff = tree_generator.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)



if __name__ == '__main__':
    unittest.main()
