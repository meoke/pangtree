import unittest
import numpy as np
from ddt import ddt, data, unpack

from context import maf_reader
from context import toolkit
from context import POAGraph
from context import Source
from context import Node
from context import cons

@ddt
class MafReaderTests(unittest.TestCase):

    @data(
        ([0.3, 0.4, 0.45, 0.6, 0.7], 1, 0.65 ,0.4)
    )
    @unpack
    def test_new_cutoff(self, compatibilities, multiplier, smallest_comp_up_to_now, expected_cutoff):
        cutoff = cons._find_cutoff_new(compatibilities, multiplier, smallest_comp_up_to_now)
        self.assertEqual(expected_cutoff, cutoff)

if __name__ == '__main__':
    unittest.main()
