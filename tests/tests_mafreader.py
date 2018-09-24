import unittest
import argparse
from ddt import ddt, data, unpack

from context import mafreader
from context import metadatareader


@ddt
class CmdargsTest(unittest.TestCase):

    def setUp(self):
        self.test1_metadata = metadatareader.read("Files/test1_metadata.json")
        self.ebola_metadata = metadatareader.read("Files/ebola_metadata.json")
        self.mycoplasma_metadata = metadatareader.read("Files/mycoplasma_metadata.json")

    @data("Files/test1.maf")
    @unittest.skip("simple")
    def test_read_maf(self, maf_path):
        graph, paths_manager = mafreader.read(maf_path, self.test1_metadata)    \

    @data("Files/ebola_ncbi.maf")
    # @unittest.skip("ebola")
    def test_read_maf_ebola(self, maf_path):
        graph, paths_manager = mafreader.read(maf_path, self.ebola_metadata)    \

    @data("Files/mycoplasma_full.maf")
    @unittest.skip("myco")
    def test_read_maf_myco(self, maf_path):
        graph, paths_manager = mafreader.read(maf_path, self.mycoplasma_metadata)


if __name__ == '__main__':
    unittest.main()
