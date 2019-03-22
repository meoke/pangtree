import unittest

from ddt import ddt, data, unpack
from tests.context import MAX1, MAX2, NODE1, NODE2, NODE3, NODE4
import numpy as np

@ddt
class NodeCutoffStrategiesTests(unittest.TestCase):

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
        node1_strategy: NODE1 = NODE1(multiplier)
        actual_cutoff = node1_strategy.find_node_cutoff(compatibilites, []).cutoff
        self.assertEqual(actual_cutoff, expected_cutoff)

    def test_node1_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node1_strategy: NODE1 = NODE1(1)
            _ = node1_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

    @data(
        (0.2, [0.1, 0.25, 0.2], [0.9], 1)
    )
    @unpack
    def test_node2_strategy_guard_greater_than_all_comps(self, expected_cutoff, compatibilities, so_far_cutoffs, multiplier):
        node2_strategy = NODE2(multiplier)
        actual_cutoff = node2_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        node1_strategy = NODE1(multiplier)
        node1_strategy_cutoff = node1_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(actual_cutoff, expected_cutoff)
        self.assertEqual(actual_cutoff, node1_strategy_cutoff)

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
        node2_strategy = NODE2(multiplier)
        actual_cutoff = node2_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff)

    def test_node2_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node2_strategy: NODE2 = NODE2(1)
            _ = node2_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")

    @data(
        # empty level guards - use max 2 strategy
        (0.7, [0.2, 0.85, 0.7, 0.8], [])
    )
    @unpack
    def test_node3_strategy_empty_level_guards(self, expected_cutoff, compatibilities, so_far_cutoffs):
        node3_strategy = NODE3()
        actual_cutoff = node3_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        max2_strategy = MAX2()
        max2_strategy_cutoff = max2_strategy.find_max_cutoff(compatibilities).cutoff
        self.assertEqual(actual_cutoff, expected_cutoff)
        self.assertEqual(actual_cutoff, max2_strategy_cutoff)

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
        node3_strategy = NODE3()
        actual_cutoff = node3_strategy.find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff)

    def test_node3_strategy_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node3_strategy: NODE3 = NODE3()
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
        node4_strategy = NODE4()
        actual_cutoff = node4_strategy.find_node_cutoff(compatibilities, []).cutoff
        max2_strategy = MAX2()
        max2_strategy_cutoff = max2_strategy.find_max_cutoff(compatibilities).cutoff
        self.assertEqual(actual_cutoff, expected_cutoff)
        self.assertEqual(actual_cutoff, max2_strategy_cutoff)

    def test_node4_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            node3_strategy: NODE4 = NODE4()
            _ = node3_strategy.find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Cannot find cutoff.")


if __name__ == '__main__':
    unittest.main()
