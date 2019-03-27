import unittest
from ddt import ddt, data, unpack

from tests.context import PangraphBuilderFromMAF
from tests.Pangraph_Tests.Pangraph_Tests import PangraphTests


@ddt
class PangraphBuilderFromMAFTest_HelpMethods(PangraphTests):

    @data(("Pangraph_Tests/PangraphBuilderFromMAF_Tests/files_help_methods/test_1_no_spaces_in_blocks.maf", 23),
          ("Pangraph_Tests/PangraphBuilderFromMAF_Tests/files_help_methods/test_2_spaces_in_blocks.maf", 18))
    @unpack
    def test_calc_nodes(self, maf_path, expected_nodes_count):
        mafalignment = [*PangraphTests.read_maf(maf_path)]
        actual_nodes_count = PangraphBuilderFromMAF.get_nodes_count(mafalignment)
        self.assertEqual(expected_nodes_count, actual_nodes_count)


if __name__ == '__main__':
    unittest.main()
