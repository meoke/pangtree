import unittest
from ddt import ddt, data, unpack

from context import tree


@ddt
class TreeTestNodeCutoff(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_node_cutoff_no_compatibilities(self):
        #moved
        with self.assertRaises(ValueError) as err:
            actual_cutoff = tree.find_node_cutoff([], 0, [0, 1])
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Finding max cutoff is not possible.")

    @data((0.5, [0.5], [], 1),
          (0.5, [0.5], [], 0),
          (0.5, [0.5], [], 0.5))
    @unpack
    def test_get_node_cutoff_single_comp_value(self, expected_cutoff, compatibilities, old_comps, multiplier):
        #moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((0.7, [0.5, 0.7], [], 1),
          (1, [1, 0.45], [], 0),
          (0.9, [0.9, 0.5], [], 0.4))
    @unpack
    def test_get_node_cutoff_two_comp_values(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((.8, [.3, .4, .8], [], 1),
          (0.91, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97], [], 1),
          (0.91, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97], [], 1),
          (0.5, [0.5], [], 1),
          (0.75, [0.1, 0.75, 0.8, 0.81, 1], [], 1),
          (0.9, [0.5, 0.9, 0.99], [], 1),
          (0.8333, [1.0, 0.9444, 0.8333, 0.0556, 0.1111], [], 1))
    @unpack
    def test_get_node_cutoff_multiplier_1(self, expected_cutoff, compatibilities, old_comps, multiplier):
        #moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((.8, [.3, .4, .8], [], 1.3),
          (0.32, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97], [], 0.01))
    @unpack
    def test_get_node_cutoff_multiplier_ok(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)

    @data((.8, [.3, .4, .8], [], 10),
          (0.91, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97], [], 10))
    @unpack
    def test_get_node_cutoff_multiplier_too_big(self, expected_cutoff, compatibilities, old_comps, multiplier):
        # moved
        actual_cutoff = tree.find_node_cutoff(compatibilities, multiplier, old_comps)
        self.assertEqual(actual_cutoff, expected_cutoff)