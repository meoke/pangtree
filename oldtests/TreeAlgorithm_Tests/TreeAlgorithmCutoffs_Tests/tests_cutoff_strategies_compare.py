import unittest

from ddt import ddt, data, unpack
from tests.context import MAX1, MAX2, NODE1, NODE2, NODE3, NODE4
import numpy as np

@ddt
class CutoffStrategiesComparisonTests(unittest.TestCase):

    @data(
        (0.5, [0.5]),
        (0.4, [0.1, 0.2, 0.3, 0.4]),
        (0.1, [0.1, 0.1, 0.1, 0.1]),
        (0.998, [0.997, 0.998, 0.999]),
        (2, [1, 2, 3, 4])
        )
    @unpack
    def test_max1_max2_strategy_should_be_equal_for_full_range(self, expected_cutoff, compatibilites):
        max1_strategy: MAX1 = MAX1([0,1])
        max1_cutoff = max1_strategy.find_max_cutoff(compatibilites).cutoff
        max2_strategy: MAX2 = MAX2()
        max2_cutoff = max2_strategy.find_max_cutoff(compatibilites).cutoff
        # self.assertEqual(max1_cutoff, max2_cutoff)
        self.assertEqual(expected_cutoff, max1_cutoff)


if __name__ == '__main__':
    unittest.main()
