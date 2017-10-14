from Multialignment import Multialignment


def convert_maf_to_po(maf_file_name,
                      merge_blocks_option,
                      visualize_option = False,
                      consensus_option = False,
                        consensus_iterative = False,
                        hbmin = 0.9,
                        min_comp = 0.1,
                      fasta_option = False):

    m = Multialignment()
    m.build_multialignment_from_maf(maf_file_name, merge_blocks_option)

    if consensus_option:
        m.generate_consensus(consensus_iterative, hbmin, min_comp)

    if visualize_option:
        m.generate_visulization()

    if fasta_option:
        m.generate_fasta_files()
