from consensus.ConsensusesTree import ConsensusesTree
from consensus.FindCutoff import MAX1, MAX2, FindMaxCutoff, FindNodeCutoff, NODE1, NODE2, NODE3, NODE4
from metadata.MultialignmentMetadata import MultialignmentMetadata
from fasta_providers.FromEntrezFastaProvider import FromEntrezFastaProvider
from fasta_providers.FromZIPFastaProvider import FromZIPSystemProvider
from pangraph.Pangraph import Pangraph
from pangraph.PangraphBuilders.PangraphBuilderFromDAG import PangraphBuilderFromDAG
from pangraph.PangraphBuilders.PangraphBuilderFromMAF import PangraphBuilderFromMAF
from pangraph.PangraphBuilders.PangraphBuilderFromPO import PangraphBuilderFromPO
from pangraph.PangraphToFilesConverters.PangraphToPO import PangraphToPO
from tools import pathtools, loggingtools
from arguments.PangenomeParameters import ConsensusAlgorithm, FastaComplementationOption, MaxCutoffOption, \
    NodeCutoffOption, PangenomeParameters
from consensus.SimplePOAConsensusGenerator import SimplePOAConsensusGenerator
from consensus.TreePOAConsensusGenerator import TreePOAConsensusGenerator
from arguments.PangenomeParameters import MultialignmentFormat
import fileformats.fasta as fileformats_fasta

global_logger = loggingtools.get_logger("")


class Pangenome:
    def __init__(self, pangenome_parameters: PangenomeParameters):
        self.params = pangenome_parameters

        self.genomes_info: MultialignmentMetadata = self._build_genomes_info()
        self.pangraph: Pangraph = Pangraph(pangenome_parameters.datatype)
        self.dagmaf = None
        self.consensuses_tree: ConsensusesTree = None
        self.missing_nucleotide_symbol = self.params.missing_nucleotide_symbol

        self._config_logging()

        global_logger.info(f'Program arguments: {str(pangenome_parameters)}')

    def run(self):
        """Creates Pangraph and runs required algorithms."""
        global_logger.info("Run Pangenome...")

        self._build_pangraph_and_update_metadata()

        if self.params.output_po:
            self._generate_po_file()

        if self.params.generate_fasta:
            self._generate_fasta_files_to_directory()

        if self.params.consensus_type != ConsensusAlgorithm.NO:
            self._generate_consensus()

    def _config_logging(self):
        if self.params.verbose:
            loggingtools.add_file_handler_to_logger(self.params.output_path, "details", "details.log",
                                                    propagate=False)
            loggingtools.add_file_handler_to_logger(self.params.output_path, "", "details.log", propagate=False)

    def _build_pangraph_and_update_metadata(self) -> None:
        if self.params.multialignment_format.value == MultialignmentFormat.MAF.value:
            if self.params.raw_maf:
                self._build_from_maf()
            else:
                self._build_from_maf_firstly_converted_to_dag()
        elif self.params.multialignment_format.value == MultialignmentFormat.PO.value:
            self._build_from_po()

    def _build_from_maf_firstly_converted_to_dag(self):
        if self.params.fasta_complementation_option is FastaComplementationOption.NO:
            fasta_source = None
        elif self.params.fasta_complementation_option is FastaComplementationOption.NCBI:
            fasta_source = FromEntrezFastaProvider(self.params.email_address, self.params.cache)
        elif self.params.fasta_complementation_option is FastaComplementationOption.LOCAL:
            fasta_source = FromZIPSystemProvider(self.params.local_fasta_dirpath)
        else:
            raise Exception("Not known fasta complementation option. "
                            "Should be of type FastaComplementationOption."
                            "Cannot build pangraph.")

        sequences_names_from_maf = PangraphBuilderFromDAG.get_sequences_names(self.params.multialignment_file_content)
        self.genomes_info.feed_with_multialignment_data(sequences_names_from_maf,
                                                        MultialignmentFormat.MAF)
        self.pangraph.build_from_maf_firstly_converted_to_dag(mafcontent=self.params.multialignment_file_content,
                                                              fasta_source=fasta_source,
                                                              genomes_info=self.genomes_info,
                                                              missing_nucleotide_symbol=self.missing_nucleotide_symbol)

    def _build_from_maf(self):
        sequences_names_from_maf = PangraphBuilderFromMAF.get_sequences_names(self.params.multialignment_file_content)
        self.genomes_info.feed_with_multialignment_data(sequences_names_from_maf,
                                                        MultialignmentFormat.MAF)
        self.pangraph.build_from_maf(self.params.multialignment_file_content,
                                     self.genomes_info)

    def _build_from_po(self):
        sequences_names_from_po = PangraphBuilderFromPO.get_sequences_names(self.params.multialignment_file_content)
        self.genomes_info.feed_with_multialignment_data(sequences_names_from_po,
                                                        MultialignmentFormat.PO)
        self.pangraph.build_from_po(pocontent=self.params.multialignment_file_content,
                                    genomes_info=self.genomes_info,
                                    missing_nucleotide_symbol=self.missing_nucleotide_symbol)

    def _generate_fasta_files_to_directory(self):
        global_logger.info("Generating fasta file with all sequences...")
        pangraph_fasta_file_content = fileformats_fasta.pangraph_to_fasta(self.pangraph)
        pangraph_fasta_file_path = pathtools.get_child_path(self.params.output_path, "sequences.fasta")
        pathtools.save_to_file(pangraph_fasta_file_content, pangraph_fasta_file_path)

        if self.consensuses_tree:
            global_logger.info("Generating fasta file with all consensues sequences...")
            consensuses_fasta_file_content = fileformats_fasta.consensuses_tree_to_fasta(self.pangraph,
                                                                                         self.consensuses_tree)
            consensuses_fasta_file_path = pathtools.get_child_path(self.params.output_path, "consensues.fasta")
            pathtools.save_to_file(consensuses_fasta_file_content, consensuses_fasta_file_path)

    def _generate_consensus(self):
        output_dir = pathtools.create_child_dir(self.params.output_path, 'consensus')

        if self.params.consensus_type == ConsensusAlgorithm.SIMPLE:
            consensus_generator = SimplePOAConsensusGenerator(
                hbmin=self.params.hbmin,
                blosum_path=self.params.blosum_file_path,
            )
            self.consensuses_tree = consensus_generator.get_consensuses_tree(
                pangraph=self.pangraph,
                output_dir=output_dir,
            )
        elif self.params.consensus_type == ConsensusAlgorithm.TREE:
            consensus_generator = TreePOAConsensusGenerator(
                max_node_strategy=self._get_max_cutoff_strategy(),
                node_cutoff_strategy=self._get_node_cutoff_strategy(),
                blosum_path=self.params.blosum_file_path,
                stop=self.params.stop,
                re_consensus=self.params.re_consensus,
                p=self.params.p
            )
            self.consensuses_tree = consensus_generator.get_consensuses_tree(
                pangraph=self.pangraph,
                output_dir=output_dir,
                log_tresholds=self.params.verbose
            )
        else:
            raise Exception("Not known consensus generation algorithm option."
                            "Should be of type ConsensusAlgorithm."
                            "Cannot generate consensuses.")

    def _get_max_cutoff_strategy(self) -> FindMaxCutoff:
        if self.params.max_cutoff_option == MaxCutoffOption.MAX1:
            return MAX1(self.params.search_range)
        elif self.params.max_cutoff_option == MaxCutoffOption.MAX2:
            return MAX2()
        else:
            raise Exception("Not known max cutoff option."
                            "Should be of type MaxCutoffOption."
                            "Cannot generate consensuses.")

    def _get_node_cutoff_strategy(self) -> FindNodeCutoff:
        if self.params.node_cutoff_option == NodeCutoffOption.NODE1:
            return NODE1(self.params.multiplier)
        elif self.params.node_cutoff_option == NodeCutoffOption.NODE2:
            return NODE2(self.params.multiplier)
        elif self.params.node_cutoff_option == NodeCutoffOption.NODE3:
            return NODE3()
        elif self.params.node_cutoff_option == NodeCutoffOption.NODE4:
            return NODE4()
        else:
            raise Exception("Not known node cutoff option."
                            "Should be of type NodeCutoffOption."
                            "Cannot generate consensuses.")

    def _build_genomes_info(self):
        return MultialignmentMetadata(self.params.metadata_file_content)

    def _generate_po_file(self) -> None:
        pangraph_to_po = PangraphToPO()
        po_file_content = pangraph_to_po.get_po_file_content_from_pangraph(self.pangraph)
        po_file_path = pathtools.get_child_path(self.params.output_path, "pangraph.po")
        pathtools.save_to_file(po_file_content, po_file_path)




