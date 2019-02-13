from io import StringIO
from pathlib import Path

from graph.FastaSource import EntrezFastaSource, FastaFileSystemSource
from metadata import reader as metadatareader
from consensus.algorithm.TreeConfig import TreeConfig
from consensus.algorithm import tree as consensustree
from consensus.algorithm import simple as consensussimple
from graph.Pangraph import Pangraph
from userio.ProgramParameters import ConsensusAlgorithm, FastaComplementationOption


class Pangenome:
    def __init__(self, metadata):
        self.genomes_info = self._build_genomes_info(metadata)
        self.dagmaf = None
        self.pangraph = Pangraph()

    def build_from_maf_firstly_converted_to_dag(self, mafcontent: StringIO, fasta_complementation_option: FastaComplementationOption, fasta_dir: Path = None):
        if fasta_complementation_option == FastaComplementationOption.No:
            fasta_source = None
        elif fasta_complementation_option == FastaComplementationOption.NCBI:
            fasta_source = EntrezFastaSource()
        elif fasta_complementation_option == FastaComplementationOption.LocalFasta:
            fasta_source = FastaFileSystemSource(fasta_dir)
        else:
            raise Exception("Not known fasta complementation option. Should be none, 'ncbi' or path to fasta directory.")

        self.genomes_info.feed_with_maf_data(mafcontent)
        self.pangraph.build_from_maf_firstly_converted_to_dag(mafcontent, fasta_source, self.genomes_info)

    def build_from_maf(self, mafcontent: StringIO):
        self.genomes_info.feed_with_maf_data(mafcontent)
        self.pangraph.build_from_maf(mafcontent, self.genomes_info)

    def generate_fasta_files(self, output_dir):
        pass

    def generate_consensus(self, output_dir, consensus_type, hbmin, r, multiplier, stop, re_consensus, anti_granular):
        if consensus_type == ConsensusAlgorithm.Simple:
            self.pangraph = consensussimple.run(output_dir, self.pangraph, hbmin, self.genomes_info)
        elif consensus_type == ConsensusAlgorithm.Tree:
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
