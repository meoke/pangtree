# from .geninfo import input_reader
import geninfo.input_reader

class Pangraph:
    def __init__(self, multialignment_file, data_file):
        self.genomes_info = self._read_genomes_info(data_file)
        self.graph = self._build_graph(multialignment_file)

    def generate_fasta_files(self, output_dir):
        pass

    def generate_consensus(self, output_dir, hbmin, mincomp, r, multiplier, stop, re_consensus):
        pass

    def generate_visualization(self, output_dir):
        pass

    def _read_genomes_info(self, data_file):
        return geninfo.input_reader.read(data_file)

    def _build_graph(self, multialignment_file):
        pass