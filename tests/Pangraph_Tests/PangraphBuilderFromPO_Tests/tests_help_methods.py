import unittest
from ddt import ddt, data, unpack

from tests.context import PangraphBuilderFromPO
from tests.Pangraph_Tests.Pangraph_Tests import PangraphTests


@ddt
class PangraphBuilderFromPOTest_HelpMethods(PangraphTests):

    def test_extract_line_value(self):
        # todo tests
        pass


if __name__ == '__main__':
    unittest.main()
