import unittest
from ddt import ddt, data, unpack

from context import tree
import numpy as np

@ddt
class TreeTestNodeCutoffAntiGranular(unittest.TestCase):

    def setUp(self):
        pass

    def test_cutoff_no_compatibilities(self):
        #moved
        with self.assertRaises(ValueError) as err:
            actual_cutoff = tree.find_node_cutoff([], 0, 0, [0, 1])
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Finding max cutoff is not possible.")

    @data((0.7, np.array([0.2, 0.7, 0.8, 0.85]), [0.75], 1))
    @unpack
    def test_cutoff_smaller_than_previous(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibility_to_node_sequences=compatibilities,
                                              multiplier=multiplier,
                                              nodeid=0,
                                              mincomps=old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((0.2, [0.2, 0.7, 0.8, 0.85], [0.3], 1))
    @unpack
    def test_cutoff_take_previous(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((0.85, [0.85, 0.7], [0.83, 0.9], 1))
    @unpack
    def test_cutoff_take_next(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((0.9, [0.2, 0.25, 0.3, 0.9], [0.9, 0.95], 1))
    @unpack
    def test_cutoff_equals_previous(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)