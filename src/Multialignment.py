# import os
# from subprocess import run
import toolkit as t
# from MAFReader import MAFReader
# from POReader import POReader
# from POAGraphVisualizer import POAGraphVisualizer
# from fasta_generators import generate_source_as_fasta_from_poagraph, generate_consensus_as_fasta_from_poagraph

class Multialignment(object):
    def __init__(self):
        self.poagraphs = None
        self.poagraphs_paths = []
        self.name = None
        self.output_dir = None


    def build_multialignment_from_maf(self, maf_file_name, merge_option):
        self.name = self._get_mutlialignement_name(maf_file_name)
        self.output_dir = self._get_output_dir(maf_file_name)


    def _get_mutlialignement_name(self, input_file_name):
        return t.get_file_name_without_extension(input_file_name)


    def _get_output_dir(self, input_file_name):
        return t.create_next_sibling_dir(input_file_name, "converted")


class NoConsensusFound(Exception):
    pass