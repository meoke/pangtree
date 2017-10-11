# import os
# from subprocess import run
import toolkit as t
import maf_reader as maf_reader
# from POReader import POReader
# from POAGraphVisualizer import POAGraphVisualizer
# from fasta_generators import generate_source_as_fasta_from_poagraph, generate_consensus_as_fasta_from_poagraph

class Multialignment(object):
    def __init__(self):
        self.name = None
        self.output_dir = None
        self.poagraphs = None


    def build_multialignment_from_maf(self, maf_file_name, merge_option):
        def _get_mutlialignment_name(input_file_name):
            return t.get_file_name_without_extension(input_file_name)

        def _get_output_dir(input_file_name):
            return t.create_next_sibling_dir(input_file_name, "converted")

        self.name = _get_mutlialignment_name(maf_file_name)
        self.output_dir = _get_output_dir(maf_file_name)
        self.poagraphs = maf_reader.parse_to_poagraphs(file_name = maf_file_name,
                                                       merge_option = merge_option,
                                                       multialignment_name = self.name,
                                                       output_dir = self.output_dir)




class NoConsensusFound(Exception):
    pass