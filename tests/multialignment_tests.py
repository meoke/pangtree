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

    @unittest.skip("To implement - Check multialignment building")
    def test_build_multialignment_from_maf(self):
        test_file = self.temp_test_file
        m = Multialignment()
        m.build_multialignment_from_maf(test_file, "all")

    @unittest.skip("To implement - Check consensus generation")
    def test_consensus_generation(self):
        test_file = self.temp_test_file
        consensus_iterative = False
        hbmin = 0.9
        min_comp = 0.1

        m = Multialignment()
        m.build_multialignment_from_maf(test_file, "all")
        m.generate_consensus(consensus_iterative, hbmin, min_comp)

    def test_find_cutoff_value(self):
        compatibilities = [0.99, 0.998, 0.7, 0.1, 1, 0.82, 0.5, 0.67, 0.99, 0.23]
        comp_range = (0.6, 0.9)
        expected_cutoff_value = 0.998

        m = Multialignment()
        actual_cutoff_value = m._find_cutoff_value(compatibilities, comp_range)

        self.assertEqual(expected_cutoff_value, actual_cutoff_value)

if __name__ == '__main__':
    unittest.main()
