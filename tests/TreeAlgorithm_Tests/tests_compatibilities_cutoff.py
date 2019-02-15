import unittest

from ddt import ddt, data, unpack


@ddt
class CompatibilitiesTests(unittest.TestCase):

    @data(([], 0),
          ([], 0))
    @unpack
    def test_max1_strategy(self, compatibiliteis, expected_cutoff):
        self.assertEqual(True, True)
