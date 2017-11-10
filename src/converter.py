from Multialignment import Multialignment


def convert_maf_to_po(file_name,
                      file_format,
                      merge_blocks_option,
                      visualize_option = False,
                      consensus_option = False,
                        hbmin = 0.9,
                        min_comp = 0.1,
                        range='[0.9,1]',
                        tresholds='[1,0.9,0.8,0.7]',
                      fasta_option = False,
                      data_type='ebola'):

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
                             tresholds=tresholds)

    if visualize_option:
        m.generate_visulization()

    if fasta_option:
        m.generate_fasta_files()

    if consensus_option or visualize_option:
        m.generate_visualization(consensus_option, visualize_option)
