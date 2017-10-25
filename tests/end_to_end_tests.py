import unittest
import time
from context import converter
from context import Multialignment


class EndToEndTest(unittest.TestCase):
    #@unittest.skip("end to end")
    def test_run_full_path(self):
        start = time.clock()
        converter.convert_maf_to_po(#file_name = 'files/ebola_100th_block/ebola_100th_block.maf',
                                    #file_name='files/entire_ebola/ebola_ncbi.maf',
                                    # maf_file_name='files/mycoplasma_maf/alignment_clean.maf',
                                    # file_name= 'files/simple/simple.maf',
                                    file_name='files/entire_ebola_po/entire_ebola.po',
                                    file_format='po',
                                    merge_blocks_option = "all",
                                    consensus_option=True,
                                    consensus_iterative = True,
                                    hbmin=0.7,
                                    min_comp=0.01,
                                    data_type='ebola'
                                    )
        end = time.clock()
        print("Running time: ", time.strftime('%H:%M:%S', time.gmtime(end-start)))
        self.assertTrue(True)

    @unittest.skip("consensus from mycoplasma po")
    def test_consensus_generation_from_po(self):
        start = time.clock()
        m = Multialignment('mycoplasma')
        m.build_multialignment_from_po(po_file_name='files/mycoplasma_po/mycoplasma.po')
        m.generate_consensus(consensus_iterative=False, hbmin=0.9, min_comp=0.1)
        m.generate_visualization(consensuses_comparison=True)
        end = time.clock()
        print("Running time: ", time.strftime('%H:%M:%S', time.gmtime(end-start)))
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
