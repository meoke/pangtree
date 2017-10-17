# import os
from subprocess import run
import toolkit as t
import maf_reader as maf_reader
import po_reader as po_reader
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
        self.poagraphs = maf_reader.parse_to_poagraphs(file_path = maf_file_name,
                                                       merge_option = merge_option,
                                                       multialignment_name = self.name,
                                                       output_dir = self.output_dir)


    def generate_consensus(self, consensus_iterative, hbmin, min_comp):
        for i, p in enumerate(self.poagraphs):
            consensus_output_dir = t.create_child_dir(p.path, "consensus")
            if consensus_iterative:
                self._run_iterative_consensus_generation(consensus_output_dir, i, hbmin, min_comp)
            else:
                self._run_single_consensus_generation(consensus_output_dir, i, hbmin)
            self._save_consensuses(consensus_output_dir, i)


    def _run_single_consensus_generation(self, consensus_output_dir, i, hbmin, consensus_name = "consensus"):
        poagraph_as_po = self.poagraphs[i].generate_po()

        file_name = t.join_path(consensus_output_dir, consensus_name)
        hb_file_name = t.change_file_extension(file_name,  '.hb')
        with open(file_name, 'w') as output_po_file:
            output_po_file.write(poagraph_as_po)
        run(['../bin/poa', '-read_msa', file_name, '-hb', '-po', hb_file_name, '../bin/blosum80.mat', '-hbmin',
             str(hbmin)])

        new_poagraph = po_reader.parse_to_poagraph(hb_file_name)

        # new_poagraph.calc_consensuses_compatibility()

    def _save_consensuses(self, output_dir, poagraph_ID):
        pass




class NoConsensusFound(Exception):
    pass