import unittest
import argparse
from ddt import ddt, data, unpack

from context import mafreader
from context import metadatareader


@ddt
class CmdargsTest(unittest.TestCase):

    def setUp(self):
        self.metadata = metadatareader.read("Files/test1_metadata.json")

    @data("Files/test1.maf")
    def test_read_maf(self, maf_path):
        graph, paths_manager = mafreader.read(maf_path, self.metadata)


if __name__ == '__main__':
    unittest.main()
