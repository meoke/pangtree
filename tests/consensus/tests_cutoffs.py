import unittest
from typing import List

from ddt import unpack, data, ddt

from tests.context import pNode, pSeq, pPoagraph, CT, P, cutoffs, Range, ConsensusInputError, Multiplier


def sid(x): return pSeq.SequenceID(x)


def nid(x): return pNode.NodeID(x)


def b(x): return pNode.Base(x)


@ddt
class CutoffsTests(unittest.TestCase):

    def setUp(self) -> None:
        self.csv_files_dir = 'tests/datamodel/input_types/csv_files/'
        self.alignment_files_dir = 'tests/datamodel/input_types/alignment_files/'

    @data((P(0.5), CT.CompatibilityToPath(0.836660026534076)),
          (P(1), CT.CompatibilityToPath(0.7)),
          (P(4), CT.CompatibilityToPath(0.6561)))
    @unpack
    def test_1_p(self, p: P, expected_cutoff: CT.CompatibilityToPath):
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
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 7, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 12, 3, 4, 5, 6, 7, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq4'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [10, 11, 2, 3, 4, 5, 6, 7, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '2'}))
        }

        pangraph = pPoagraph.Poagraph(nodes, sequences)

        consensus_path = pSeq.SequencePath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 18, 19])])
        compatibilities = pangraph.get_compatibilities(pangraph.get_sequences_ids(), consensus_path, p)
        max2_strategy: cutoffs.MAX2 = cutoffs.MAX2()
        actual_max2_cutoff = max2_strategy.find_max_cutoff([*compatibilities.values()]).cutoff
        self.assertAlmostEqual(expected_cutoff.value, actual_max2_cutoff.value)

    @data(
        (0.5, [CT.CompatibilityToPath(0.5)]),
        (0.4, [*map(CT.CompatibilityToPath, [0.1, 0.2, 0.3, 0.4])]),
        (0.1, [*map(CT.CompatibilityToPath, [0.1, 0.1, 0.1, 0.1])]),
        (0.998, [*map(CT.CompatibilityToPath, [0.997, 0.998, 0.999])]),
        (2, [*map(CT.CompatibilityToPath, [1, 2, 3, 4])])
        )
    @unpack
    def test_2_max1_max2_strategy_should_be_equal_for_full_range(self, expected_cutoff_value, compatibilites):
        max1_strategy: cutoffs.MAX1 = cutoffs.MAX1(Range([0,1]))
        max1_cutoff = max1_strategy.find_max_cutoff(compatibilites).cutoff
        max2_strategy: cutoffs.MAX2 = cutoffs.MAX2()
        max2_cutoff = max2_strategy.find_max_cutoff(compatibilites).cutoff
        self.assertEqual(max1_cutoff, max2_cutoff)
        self.assertEqual(expected_cutoff_value, max1_cutoff.value)

    @data(
        # single compatibility value
        (0.5, [0.5], [0, 1]),
        (0.5, [0.5], [0, 0]),
        (0.5, [0.5], [.1, .7]),
        (0.5, [0.5], [1, 1]),

        # two compatibilities values
        (0.7, [0.5, 0.7], [0, 1]),
        (1, [1, 0.45], [.1, .9]),
        (0.9, [0.9, 0.5], [.01, .8]),

        # repeated values
        (0.7, [0.5, 0.7, 0.7], [0, 1]),
        (0.9, [0.9, 0.5, 0.5], [0, 1]),
        (1, [1, 0.45, 0.45, 0.45], [0, 1]),

        # correct search range
        (.8, [.3, .4, .8], [0, 1]),
        (0.91, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97], [0, 0.5]),
        (0.91, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97], [0, 0.5]),
        (0.5, [0.5], [.1, .7]),
        (1, [0.1, 0.75, 0.8, 0.81, 1], [.6, 1]),
        (0.99, [0.5, 0.9, 0.99], [1, 1]),

        # single value in search range
        (.4, [.3, .4, .8], [0.4, 0.45])
        )
    @unpack
    def test_3_max1_strategy_p_1(self, expected_cutoff, compatibiliteis, cutoff_search_range):
        max1_strategy = cutoffs.MAX1(Range(cutoff_search_range))
        actual_cutoff = max1_strategy.find_max_cutoff([CT.CompatibilityToPath(c) for c in compatibiliteis]).cutoff
        self.assertEqual(actual_cutoff.value, expected_cutoff)

    def test_4_max1_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            max1_strategy = cutoffs.MAX1(Range([0, 1]))
            _ = max1_strategy.find_max_cutoff([]).cutoff
        self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

    def test_5_max_1_incorrect_search_range_order(self):
        with self.assertRaises(ConsensusInputError) as err:
            max1_strategy = cutoffs.MAX1(Range([1, 0]))
            _ = max1_strategy.find_max_cutoff([CT.CompatibilityToPath(0.2), CT.CompatibilityToPath(0.3)]).cutoff
        self.assertEqual(str(err.exception), "CUTOFF SEARCH RANGE first value must be"
                                             " smaller or equal to second value.")

    def test_6_max_1_incorrect_search_range_length(self):
        with self.assertRaises(ConsensusInputError) as err:
            max1_strategy = cutoffs.MAX1(Range([0, 0.5, 1]))
            _ = max1_strategy.find_max_cutoff([CT.CompatibilityToPath(0.2), CT.CompatibilityToPath(0.3)]).cutoff
        self.assertEqual(str(err.exception), "CUTOFF SEARCH RANGE must have length 2.")

    @data(
        # single compatibility value
        (0.5, [CT.CompatibilityToPath(0.5)]),
        (0.5, [CT.CompatibilityToPath(0.5)]),
        (0.5, [CT.CompatibilityToPath(0.5)]),
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

        # repeated distance between values
        (.4, [*map(CT.CompatibilityToPath, [.3, .4, .5])]),

        # all the same values
        (.1, [*map(CT.CompatibilityToPath, [.1, .1, .1])])
    )
    @unpack
    def test_max2_strategy(self, expected_cutoff: float, compatibilities: List[CT.CompatibilityToPath]):
        max2_strategy = cutoffs.MAX2()
        actual_cutoff = max2_strategy.find_max_cutoff(compatibilities).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)

    @data(
        # single compatibility value
        (0.5, [0.5], 1),
        (0.5, [0.5], 0),
        (0.5, [0.5], 0.5),

        # two compatibilities values
        (0.7, [0.5, 0.7], 1),
        (1, [1, 0.45], 0),
        (0.9, [0.9, 0.5], 0.4),

        # multiplier equals 1
        (.8, [.3, .4, .8], 1),
        (0.91, [0.31, 0.97, 0.32, 0.91, 0.92, 0.93], 1),
        (0.91, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97], 1),
        (0.5, [0.5], 1),
        (0.75, [0.1, 1, 0.75, 0.8, 0.81], 1),
        (0.9, [0.5, 0.9, 0.99], 1),
        (0.8333, [1.0, 0.9444, 0.8333, 0.0556, 0.1111], 1),

        # correct multiplier, not 1
        (.8, [.3, .4, .8], 1.3),
        (0.32, [0.31, 0.97, 0.32, 0.91, 0.92, 0.93], 0.01),

        # incorrect multiplier, too big
        (.8, [.3, .4, .8], 10),
        (0.91, [0.31, 0.97, .32, 0.91, 0.92, 0.93], 10)
        )
    @unpack
    def test_node1_strategy(self, expected_cutoff, compatibilites, multiplier):
        node1_strategy: cutoffs.NODE1 = cutoffs.NODE1(Multiplier(multiplier))
        compatibilites = [CT.CompatibilityToPath(c) for c in compatibilites]
        actual_cutoff = node1_strategy.find_node_cutoff(compatibilites, []).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)

    def test_node1_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node1_strategy: cutoffs.NODE1 = cutoffs.NODE1(Multiplier(1))
            _ = node1_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

    @data(
        (0.2, [0.1, 0.25, 0.2], [0.9], 1)
    )
    @unpack
    def test_node2_strategy_guard_greater_than_all_comps(self, expected_cutoff, compatibilities, so_far_cutoffs, multiplier):
        node2_strategy = cutoffs.NODE2(Multiplier(multiplier))
        compatibilities = [CT.CompatibilityToPath(c) for c in compatibilities]
        so_far_cutoffs = [CT.CompatibilityToPath(c) for c in so_far_cutoffs]
        actual_cutoff = node2_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        node1_strategy = cutoffs.NODE1(Multiplier(multiplier))
        node1_strategy_cutoff = node1_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)
        self.assertEqual(node1_strategy_cutoff.value, actual_cutoff.value)

    @data(
        # guard < all compatibilities
        #
        # take smaller then so far cutoffs
        (0.7, [0.2, 0.7, 0.8, 0.85], [0.75], 1),
        (0.2, [0.2, 0.85, 0.7, 0.8], [0.3], 1),

        # take smaller then so far cutoffs
        (0.7, [0.85, 0.7], [0.83, 0.9], 1),

        # take equals to min of so far cutoffs
        (0.9, [0.2, 0.9, 0.25, 0.3], [0.9, 0.95], 1)
    )
    @unpack
    def test_node2_strategy(self, expected_cutoff, compatibilities, so_far_cutoffs, multiplier):
        node2_strategy = cutoffs.NODE2(Multiplier(multiplier))
        compatibilities = [CT.CompatibilityToPath(c) for c in compatibilities]
        so_far_cutoffs = [CT.CompatibilityToPath(c) for c in so_far_cutoffs]
        actual_cutoff = node2_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)

    def test_node2_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node2_strategy: cutoffs.NODE2 = cutoffs.NODE2(1)
            _ = node2_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

    @data(
        # empty level guards - use max 2 strategy
        (0.7, [0.2, 0.85, 0.7, 0.8], [])
    )
    @unpack
    def test_node3_strategy_empty_level_guards(self, expected_cutoff, compatibilities, so_far_cutoffs):
        node3_strategy = cutoffs.NODE3()
        compatibilities = [CT.CompatibilityToPath(c) for c in compatibilities]
        so_far_cutoffs = [CT.CompatibilityToPath(c) for c in so_far_cutoffs]
        actual_cutoff = node3_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        max2_strategy = cutoffs.MAX2()
        max2_strategy_cutoff = max2_strategy.find_max_cutoff(compatibilities).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)
        self.assertEqual(max2_strategy_cutoff.value, actual_cutoff.value)

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
    def test_node3_strategy(self, expected_cutoff, compatibilities, so_far_cutoffs):
        node3_strategy = cutoffs.NODE3()
        compatibilities = [CT.CompatibilityToPath(c) for c in compatibilities]
        so_far_cutoffs = [CT.CompatibilityToPath(c) for c in so_far_cutoffs]
        actual_cutoff = node3_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)

    def test_node3_strategy_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node3_strategy: cutoffs.NODE3 = cutoffs.NODE3()
            _ = node3_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

    @data(
        # empty level guard
        (0.7, [0.7, 0.8, 0.85, 0.2]),
        (0.99, [0.99, 0.9, 0.99]),
        (0.99, [0.99])
    )
    @unpack
    def test_node4_strategy_empty_level_guards(self, expected_cutoff, compatibilities):
        compatibilities = [CT.CompatibilityToPath(c) for c in compatibilities]
        node4_strategy = cutoffs.NODE4()
        actual_cutoff = node4_strategy.find_node_cutoff(compatibilities, []).cutoff
        max2_strategy = cutoffs.MAX2()
        max2_strategy_cutoff = max2_strategy.find_max_cutoff(compatibilities).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)
        self.assertEqual(max2_strategy_cutoff.value, actual_cutoff.value)

    def test_node4_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node3_strategy: cutoffs.NODE4 = cutoffs.NODE4()
            _ = node3_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

if __name__ == '__main__':
    unittest.main()