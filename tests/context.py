import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangenome')))


from pangenome.datamodel import Node as pNode
from pangenome.datamodel import Poagraph as pPoagraph
from pangenome.datamodel import Sequence as pSeq
from pangenome.datamodel import Poagraph as pPoagraph
from pangenome.datamodel.builders import maf2poagraph, dagmaf2poagraph
from pangenome.datamodel.input_types import Maf, MetadataCSV
import pangenome.tools.path as pathtools
from pangenome.datamodel.fasta_providers.ConstSymbol import ConstSymbol
from pangenome.datamodel.fasta_providers.FastaProvider import FastaProvider
from pangenome.datamodel.input_types import MissingSymbol, InputError

# from pangenome.pang.arguments.PangenomeParameters import MultialignmentFormat
#
# from metadata import MultialignmentMetadata
#
# from pangenome.pang.pangraph.Node import Node
# from pangenome.pang.pangraph.Pangraph import Pangraph
# from pangenome.pang.fasta_providers.FastaProvider import FastaProvider
# from pangenome.pang.fasta_providers.FromFileFastaProvider import FromFileFastaProvider
# from pangenome.pang.fasta_providers.FromEntrezFastaProvider import FromEntrezFastaProvider
# from pangenome.pang.pangraph.PangraphBuilders.PangraphBuilderFromDAG import PangraphBuilderFromDAG
# from pangenome.pang.pangraph.PangraphBuilders.PangraphBuilderFromMAF import PangraphBuilderFromMAF
# from pangenome.pang.pangraph.custom_types import SequenceID
#
# from pangenome.pang.fileformats.maf.reader import maf_to_dagmaf
# import pangenome.pang.fileformats.fasta as fasta
#
# from pangenome.pang.consensus.TreePOAConsensusGenerator import TreePOAConsensusGenerator
# from pangenome.pang.consensus.top_consensus import get_top_consensus, PangraphPOTranslator
#
# from pangenome.pang.pangraph.CompatibilityToPath import CompatibilityToPath
# from pangenome.pang.pangraph.DataType import DataType
# from pangenome.pang.consensus.ConsensusNode import ConsensusNode, ConsensusNodeID
# from pangenome.pang.consensus.ConsensusesTree import ConsensusesTree
#
# from pangenome.pang.pangraph.Node import NodeID
# from pangenome.pang.tools import pathtools
#
