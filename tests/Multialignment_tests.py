import unittest

from context import Multialignment
from context import toolkit


class MultialignmentTests(unittest.TestCase):
    def setUp(self):
        test_file_path = 'files/ebola_100th_block/ebola_100th_block.maf'
        self.temp_dir = toolkit.create_next_sibling_dir(test_file_path, 'test')
        self.temp_test_file = toolkit.copy_file_to_dir(file_path = test_file_path,
                                 destination_dir = self.temp_dir)

    def tearDown(self):
        toolkit.remove_dir(self.temp_dir)


    def test_build_multialignment_from_maf(self):
        test_file = self.temp_test_file
        m = Multialignment()
        m.build_multialignment_from_maf(test_file, "all")
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
