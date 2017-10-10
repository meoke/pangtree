import unittest
from context import converter

class EndToEndTest(unittest.TestCase):
    def test_run_full_path(self):
        converter.convert_maf_to_po(maf_file_name = 'temp/entire_ebola/ebola_ncbi.maf',
                                    merge_blocks_option = "all")

        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
