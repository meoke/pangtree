import unittest
from context import converter

class EndToEndTest(unittest.TestCase):
    def test_run_full_path(self):
        converter.convert_maf_to_po(maf_file_name = 'files/ebola_100th_block/ebola_100th_block.maf',
                                    merge_blocks_option = "all")

        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
