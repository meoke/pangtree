import unittest
from context import converter
from context import Multialignment


class EndToEndTest(unittest.TestCase):
    #@unittest.skip("end to end")
    def test_run_full_path(self):

        converter.convert_maf_to_po(file_name = 'files/ebola_100th_block/ebola_100th_block.maf', #it's a short file, good for testing
                                    #file_name='files/entire_ebola/ebola_ncbi.maf',
                                    #file_name='files/mycoplasma_maf/alignment_clean.maf',
                                    #file_name='files/mycoplasma_maf/alignment.maf',
                                    #file_name='files/mycoplasma2/alignment.maf',
                                    # file_name='files/mycoplasma2/alignment_bez_GCA_002205575.maf',
                                    #file_name= 'files/simple/simple.maf',
                                    #file_name='files/entire_ebola_po/entire_ebola.po',
                                    file_format='maf',
                                    merge_blocks_option="all",
                                    draw_poagraph_option=True,
                                    consensus_option=3,
                                        #min_comp=0.01,
                                        range='[90,100]',
                                        multiplier=1,
                                        stop=0.8,
                                        re_consensus = True,
                                    fasta_option=False,
                                    data_type='ebola',
                                    blocks_option=True
                                    )

        self.assertTrue(True)


    @unittest.skip("consensus from mycoplasma po")
    def test_consensus_generation_from_po(self):
        m = Multialignment('mycoplasma')
        m.build_multialignment_from_po(po_file_name='files/mycoplasma_po/mycoplasma.po')
        m.generate_consensus(consensus_iterative=False, hbmin=0.9, min_comp=0.1)
        m.generate_visualization(consensuses_comparison=True)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
