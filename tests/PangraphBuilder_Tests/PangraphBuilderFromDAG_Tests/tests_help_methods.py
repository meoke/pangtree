import unittest
from ddt import ddt, data, unpack

from context import maf_to_dagmaf

from graph.Pangraph import PangraphBuilderFromDAG
from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests


@ddt
class PangraphBuilderFromDAGTest_HelpMethods(PangraphBuilderTests):

    @data(("PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_1_missing_sequence_start.maf", 12))
    @unpack
    def test_calc_nodes(self, maf_path, expected_nodes_count):
        dagmaf = maf_to_dagmaf(maf_path)
        actual_nodes_count = PangraphBuilderFromDAG.get_nodes_count(dagmaf)
        self.assertEqual(expected_nodes_count, actual_nodes_count)

    @data(("PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_help_methods/test_1_parallel_blocks_1st_and_2nd_merge_into_3rd.maf", [0, 1]),
          ("PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_help_methods/test_2_two_blocks_in_disconnected_graph.maf", [0, 1]),
          ("PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_0_simple.maf", [0]))
    @unpack
    def test_get_starting_blocks(self, maf_path, expected_blocks_ids):
        dagmaf = maf_to_dagmaf(maf_path)
        actual_starting_blocks = PangraphBuilderFromDAG.get_starting_blocks(dagmaf)
        self.assertEqual(expected_blocks_ids, actual_starting_blocks)


if __name__ == '__main__':
    unittest.main()
