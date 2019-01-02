import metadata.reader as metadatareader
from graph import mafreader
# from consensus import algorithm as consensussimple, algorithm as consensustree
from consensus.algorithm.TreeConfig import TreeConfig
from consensus.algorithm import tree as consensustree
from consensus.algorithm import simple as consensussimple


class Pangenome:
    def __init__(self, multialignment, metadata):
        self.genomes_info = self._build_genomes_info(metadata)
        self.pangraph = self._build_pangraph(multialignment)

    def generate_fasta_files(self, output_dir):
        pass

    def generate_consensus(self, output_dir, consensus_type, hbmin, r, multiplier, stop, re_consensus, anti_granular):
        if consensus_type == 'simple':
            self.pangraph = consensussimple.run(output_dir, self.pangraph, hbmin, self.genomes_info)
        elif consensus_type == 'tree':
            tree_config = TreeConfig(hbmin=hbmin,
                                     r=r,
                                     multiplier=multiplier,
                                     stop=stop,
                                     re_consensus=re_consensus,
                                     anti_granular=anti_granular)
            self.pangraph = consensustree.run(output_dir, self.pangraph, tree_config, self.genomes_info)

    # def generate_visualization(self, output_dir):
    #     pass

    def _build_genomes_info(self, genomes_metadata):
        #todo check if metadata given else generate
        return metadatareader.read(genomes_metadata)

    def _build_pangraph(self, multialignment):
        return mafreader.read(multialignment, self.genomes_info)
