import unittest
import argparse
from ddt import ddt, data, unpack

from context import tree


@ddt
class ConsensusTreeTest(unittest.TestCase):

    @data(([0, 0.25, 0.5, 0.76, 1], [0, 100], 0.76),
          ([0, 0.25, 0.5, 0.76, 1], [75, 100], 1))
    @unpack
    def test_find_max_cutoff(self, compatibilities, cutoff_search_range, expected_cutoff):
        actual_cutoff = tree.find_max_cutoff(compatibilities, cutoff_search_range)

        self.assertEqual(actual_cutoff, expected_cutoff)



if __name__ == '__main__':
    unittest.main()