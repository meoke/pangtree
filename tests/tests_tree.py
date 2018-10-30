import unittest
from ddt import ddt, data, unpack

from context import tree


@ddt
class TreeTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_get_max_cutoff_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            actual_cutoff = tree.find_max_cutoff([], [0, 1])
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Finding max cutoff is not possible.")

    def test_get_max_cutoff_incorrect_search_range_order(self):
        with self.assertRaises(ValueError) as err:
            actual_cutoff = tree.find_max_cutoff([0.2, 0.3], [1, 0])
            self.assertEqual(str(err.exception), "For cutoff search range [x, y] x must be <= y.")

    def test_get_max_cutoff_incorrect_search_range_length(self):
        with self.assertRaises(ValueError) as err:
            actual_cutoff = tree.find_max_cutoff([0.2, 0.3], [0, 0.5, 1])
            self.assertEqual(str(err.exception), "Cutoff search range must have length 2.")

    @data((0.5, [0.5], [0, 1]),
          (0.5, [0.5], [0, 0]),
          (0.5, [0.5], [.1, .7]),
          (0.5, [0.5], [1, 1]))
    @unpack
    def test_get_max_cutoff_single_comp_value(self, expected_cutoff, compatibilities, cutoff_search_range):

        actual_cutoff = tree.find_max_cutoff(compatibilities, cutoff_search_range)
        self.assertEquals(actual_cutoff, expected_cutoff)

    @data((0.5, [0.5, 0.7], [0, 1]),
          (0.45, [1, 0.45], [.1, .9]),
          (0.5, [0.9, 0.5], [.01, .8]))
    @unpack
    def test_get_max_cutoff_two_comp_values(self, expected_cutoff, compatibilities, cutoff_search_range):
        actual_cutoff = tree.find_max_cutoff(compatibilities, cutoff_search_range)
        self.assertEquals(actual_cutoff, expected_cutoff)

    @data((.8, [.3, .4, .8], [0, 1]),
          (0.91, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97], [0, 0.5]),
          (0.91, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97], [0, 0.5]),
          (0.5, [0.5], [.1, .7]),
          (1, [0.1, 0.75, 0.8, 0.81, 1], [.6, 1]),
          (0.99, [0.5, 0.9, 0.99], [1, 1]))
    @unpack
    def test_get_max_cutoff_correct_search_range(self, expected_cutoff, compatibilities, cutoff_search_range):
        actual_cutoff = tree.find_max_cutoff(compatibilities, cutoff_search_range)
        self.assertEquals(actual_cutoff, expected_cutoff)

    @data((.4, [.3, .4, .8], [0.4, 0.45]))
    @unpack
    def test_get_max_cutoff_single_value_in_search_range(self, expected_cutoff, compatibilities, cutoff_search_range):
        actual_cutoff = tree.find_max_cutoff(compatibilities, cutoff_search_range)
        self.assertEquals(actual_cutoff, expected_cutoff)

    # @data((, []),
    #       (, [0.5]),
    #       (, [0, 2, 3]))
    # @unpack
    # def test_get_node_cutoff(self, expected_cutoff, compatibilities):
    #     actual_cutoff = tree.find_node_cutoff(compatibilities)
    #     self.assertEquals(actual_cutoff, expected_cutoff)