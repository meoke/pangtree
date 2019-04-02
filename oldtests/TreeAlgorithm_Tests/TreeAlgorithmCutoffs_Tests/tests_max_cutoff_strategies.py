import unittest

from ddt import ddt, data, unpack
from tests.context import MAX1, MAX2, NODE1, NODE2, NODE3, NODE4


@ddt
class MaxCutoffStrategiesTests(unittest.TestCase):

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
    def test_max1_strategy_p_1(self, expected_cutoff, compatibiliteis, cutoff_search_range):
        max1_strategy = MAX1(cutoff_search_range)
        actual_cutoff = max1_strategy.find_max_cutoff(compatibiliteis).cutoff
        self.assertEqual(actual_cutoff, expected_cutoff)

    def test_max1_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            max1_strategy = MAX1([0, 1])
            _ = max1_strategy.find_max_cutoff([]).cutoff
            self.assertEqual(str(err.exception), f"Empty compatibilities list. Finding max cutoff is not possible.")

    def test_max_1_incorrect_search_range_order(self):
        with self.assertRaises(ValueError) as err:
            max1_strategy = MAX1([1, 0])
            _ = max1_strategy.find_max_cutoff([0.2, 0.3]).cutoff
            self.assertEqual(str(err.exception), "For cutoff search range [x, y] x must be <= y.")

    def test_max_1_incorrect_search_range_length(self):
        with self.assertRaises(ValueError) as err:
            max1_strategy = MAX1([0, 0.5, 1])
            _ = max1_strategy.find_max_cutoff([0.2, 0.3]).cutoff
            self.assertEqual(str(err.exception), "Cutoff search range must have length 2.")

    @data(
        # single compatibility value
        (0.5, [0.5]),
        (0.5, [0.5]),
        (0.5, [0.5]),
        (0.5, [0.5]),

        # two compatibilities values
        (0.7, [0.5, 0.7]),
        (1, [1, 0.45]),
        (0.9, [0.9, 0.5]),

        # repeated values
        (0.7, [0.5, 0.7, 0.7]),
        (0.9, [0.9, 0.5, 0.5]),
        (1, [0.45, 1, 0.45, 0.45]),

        # many unique compatibilities values
        (.8, [.3, .4, .8]),
        (0.91, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97]),
        (0.91, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97]),
        (1, [0.81, 0.75, 0.8, 0.81, 1]),
        (0.9, [0.5, 0.9, 0.99]),

        # repeated distance between values
        (.4, [.3, .4, .5]),

        # all the same values
        (.1, [.1, .1, .1])
    )
    @unpack
    def test_max2_strategy(self, expected_cutoff, compatibilities):
        max2_strategy = MAX2()
        actual_cutoff = max2_strategy.find_max_cutoff(compatibilities).cutoff
        self.assertEqual(actual_cutoff, expected_cutoff)


if __name__ == '__main__':
    unittest.main()
