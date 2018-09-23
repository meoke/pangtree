import metadata.reader
from .graph import mafreader


class Pangraph:
    def __init__(self, multialignment_file, data_file):
        self.genomes_info = self._read_genomes_info(data_file)
        self.graph = None
        self.paths_manager = None
        self._build_graph(multialignment_file)

    def generate_fasta_files(self, output_dir):
        pass

    def generate_consensus(self, output_dir, hbmin, mincomp, r, multiplier, stop, re_consensus):
        # to będą info o przechodzących przez graf sekwencjach i
        pass

    def generate_visualization(self, output_dir):
        pass

    def _read_genomes_info(self, data_file):
        return metadata.reader.read(data_file)

    def _build_graph(self, multialignment_file):
        self.graph, self.paths_manager = mafreader.read(multialignment_file)
