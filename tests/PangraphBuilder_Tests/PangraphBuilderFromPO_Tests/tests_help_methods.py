import unittest
from ddt import ddt, data, unpack

from tests.context import PangraphBuilderFromPO
from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests


@ddt
class PangraphBuilderFromPOTest_HelpMethods(PangraphBuilderTests):

    def test_extract_line_value(self):
        pass


if __name__ == '__main__':
    unittest.main()
