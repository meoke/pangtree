import metadata.reader
from .graph import mafreader
import consensus.simple
import consensus.tree

class Pangenome:
    def __init__(self, multialignment_file, data_file):
        self.genomes_info = self._read_genomes_info(data_file)
        self.pangraph = None
        self._build_graph(multialignment_file)

    def generate_fasta_files(self, output_dir):
        pass

    def generate_consensus(self, output_dir, consensus_type, hbmin, mincomp, r, multiplier, stop, re_consensus):
        if consensus_type == 'simple':
            self.pangraph = consensus.simple.run(output_dir, self.pangraph, hbmin, self.genomes_info)
        else:
            self.pangraph = consensus.tree.run(self.pangraph)

    def generate_visualization(self, output_dir):
        pass

    def _read_genomes_info(self, data_file):
        return metadata.reader.read(data_file)

    def _build_graph(self, multialignment_file):
        self.pangraph = mafreader.read(multialignment_file, self.genomes_info)
