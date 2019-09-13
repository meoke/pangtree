import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))


from pangtreebuild.datamodel import Node as pNode
from pangtreebuild.datamodel import Poagraph as pPoagraph
from pangtreebuild.datamodel import Sequence as pSeq
from pangtreebuild.datamodel.Sequence import SequenceID
from pangtreebuild.datamodel import Poagraph as pPoagraph
from pangtreebuild.datamodel.builders import maf2poagraph, dagmaf2poagraph, po2poagraph
from pangtreebuild.datamodel.input_types import Maf, MetadataCSV
import pangtreebuild.tools.pathtools as pathtools
from pangtreebuild.datamodel.fasta_providers.ConstSymbolProvider import ConstSymbolProvider
from pangtreebuild.datamodel.fasta_providers.FromNCBI import FromNCBI

from pangtreebuild.datamodel.fasta_providers.FastaProvider import FastaProvider, FastaProviderException
from pangtreebuild.datamodel.fasta_providers.FromFile import FromFile
from pangtreebuild.datamodel.input_types import MissingSymbol, InputError, Po

from pangtreebuild.output import PangenomePO

from pangtreebuild.output import PangenomeFASTA
from pangtreebuild.consensus import ConsensusTree as CT
from pangtreebuild.consensus.input_types import P, Range, ConsensusInputError, Multiplier
from pangtreebuild.consensus import cutoffs, poa


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
