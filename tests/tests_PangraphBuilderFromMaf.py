import unittest
from ddt import ddt, data, unpack

from context import maf_to_dagmaf

from graph.Pangraph import PangraphBuilderFromDAG


@ddt
class PangraphBuilderFromMafTest(unittest.TestCase):

    @data(("Files/maf_gaps/test_1_left.maf", 12))
    @unpack
    def test_calc_nodes(self, maf_path, expected_nodes_count):
        dagmaf = maf_to_dagmaf(maf_path)
        actual_nodes_count = PangraphBuilderFromDAG.get_nodes_count(dagmaf)
        self.assertEqual(expected_nodes_count, actual_nodes_count)


if __name__ == '__main__':
    unittest.main()
