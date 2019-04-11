import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangenome')))


from pangpang.datamodel import Node as pNode
from pangpang.datamodel import Poagraph as pPoagraph
from pangpang.datamodel import Sequence as pSeq
from pangpang.datamodel import Poagraph as pPoagraph
from pangpang.datamodel.builders import maf2poagraph, dagmaf2poagraph, po2poagraph
from pangpang.datamodel.input_types import Maf, MetadataCSV
import pangpang.tools.pathtools as pathtools
from pangpang.datamodel.fasta_providers.ConstSymbolProvider import ConstSymbolProvider
from pangpang.datamodel.fasta_providers.FromNCBI import FromNCBI, EmailAddress

from pangpang.datamodel.fasta_providers.FastaProvider import FastaProvider, FastaProviderException
from pangpang.datamodel.fasta_providers.FromFile import FromFile
from pangpang.datamodel.input_types import MissingSymbol, InputError, Po

from pangpang.output import PangenomePO

from pangpang.output import PangenomeFASTA
from pangpang.consensus import ConsensusTree as CT
from pangpang.consensus.input_types import P


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
