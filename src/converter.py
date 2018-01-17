import time
from Multialignment import Multialignment


def convert_maf_to_po(file_name,
                      file_format,
                      merge_blocks_option,
                      draw_poagraph_option=False,
                      consensus_option=False,
                        hbmin=0.9,
                        min_comp=0.1,
                        range='[0.9,1]',
                        multiplier=1,
                        stop=0.99,
                        re_consensus = True,
                        tresholds='[1,0.9,0.8,0.7]',
                      fasta_option=False,
                      data_type='ebola'):
    start = time.clock()
    m = Multialignment(data_type)
    if file_format == 'maf':
        m.build_multialignment_from_maf(file_name, merge_blocks_option)
    elif file_format == 'po':
        m.build_multialignment_from_po(file_name)

    if consensus_option:
        m.generate_consensus(option=consensus_option,
                             hbmin=hbmin,
                             min_comp=min_comp,
                             comp_range=range,
                             tresholds=tresholds,
                             multiplier=multiplier,
                             stop=stop,
                             re_consensus=re_consensus)
    if fasta_option:
        m.generate_fasta_files()

    end = time.clock()
    processing_time = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print(processing_time)
    if consensus_option or draw_poagraph_option:
        m.generate_visualization(consensus_option, draw_poagraph_option, processing_time, m.tresholds, consensus_option)




