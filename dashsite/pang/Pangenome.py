from io import StringIO

from metadata import reader as metadatareader
from graph import mafreader
from consensus.algorithm.TreeConfig import TreeConfig
from consensus.algorithm import tree as consensustree
from consensus.algorithm import simple as consensussimple
from fileformats.maf.reader import maf_to_dagmaf
from graph.Pangraph import Pangraph


class Pangenome:
    def __init__(self, metadata):
        self.genomes_info = self._build_genomes_info(metadata)
        self.pangraph = Pangraph()

    def build_from_maf_firstly_converted_to_dag(self, mafcontent: StringIO, fasta_complementation_option):
        dagmaf = maf_to_dagmaf(mafcontent)
        self.genomes_info.feed_with_maf_data(mafcontent)
        self.pangraph.build_from_dag(dagmaf, fasta_complementation_option, self.genomes_info)

    def build_from_maf(self, mafcontent: StringIO):
        self.genomes_info.feed_with_maf_data(mafcontent)
        self.pangraph.build_from_maf(mafcontent, self.genomes_info)

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

    # def _build_pangraph(self, multialignment):
    #     #todo it should be built based on blocks not raw multialignment
    #     return mafreader.read(multialignment, self.genomes_info)

    # def _build_pangraph(self):
    #     #todo it should be built based on blocks not raw multialignment
    #     # return mafreader.read(self.dagmaf, self.genomes_info)
    #     return dagmaf_to_pangraph(self.dagmaf, self.genomes_info)
