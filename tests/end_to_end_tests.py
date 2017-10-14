import unittest
import time
from context import converter


class EndToEndTest(unittest.TestCase):
    def test_run_full_path(self):
        start = time.clock()
        converter.convert_maf_to_po(maf_file_name = 'files/ebola_100th_block/ebola_100th_block.maf',
                                    #maf_file_name='files/entire_ebola/ebola_ncbi.maf',
                                    merge_blocks_option = "all",
                                    consensus_option=True,
                                    consensus_iterative = False,
                                    hbmin=0.9,
                                    min_comp=0.1,
                                    )
        end = time.clock()
        print("Running time: ", str(end-start))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
