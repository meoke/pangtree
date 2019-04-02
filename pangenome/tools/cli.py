import argparse
import inspect
from io import StringIO
from pathlib import Path
from typing import TypeVar, Callable, Optional

from datamodel.DataType import DataType
from datamodel.builders import PoagraphBuildException
from datamodel.fasta_providers.FastaProvider import FastaProviderOption, EmailAddress
from datamodel.input_types import Maf, MetadataCSV


class InvalidPath(Exception):
    pass


def _get_path_if_valid(path: str) -> Path:
    """Check if path exists."""

    file_path = Path(path)
    if not file_path.is_file():
        raise InvalidPath(f"{file_path}")
    return file_path


def _data_type(data_type: str) -> DataType:
    """Converts command line argument to DataType"""
    if data_type == "p":
        data_type_full_name = "Proteins"
    elif data_type == "n":
        data_type_full_name = "Nucleotides"
    else:
        raise argparse.ArgumentError("Unknown data type. \'p\' for proteins or \'n\' for nucleotides available.")
    try:
        dt = DataType[data_type_full_name]
        return dt
    except KeyError:
        raise argparse.ArgumentError("Data type parsing error.")


def _fasta_provider_option(arg_fasta_provider_option: str) -> FastaProviderOption:
    """Converts command line argument to FastaProviderOption"""

    try:
        return FastaProviderOption[arg_fasta_provider_option.upper()]
    except KeyError:
        raise argparse.ArgumentError("Incorrect FASTA_PROVIDER argument.")


T = TypeVar('T')


def _cli_file_arg(arg: str, constructor: Callable[[StringIO, Optional[Path]], T]) -> T:
    try:
        filepath = _get_path_if_valid(arg)
    except InvalidPath:
        raise argparse.ArgumentTypeError(f"File {arg} does not exist or is not a file.")
    with open(arg) as infile:
        file_content = StringIO(infile.read())
        try:
            return constructor(file_content, filepath)
        except PoagraphBuildException as p:
            raise argparse.ArgumentError("Incorrect file content") from p


def _maf_file(x: str) -> Maf:
    return _cli_file_arg(x, Maf)


def _metadata_file(x: str) -> MetadataCSV:
    return _cli_file_arg(x, MetadataCSV)


def cli_arg(constructor: Callable[[str], T]) -> Callable[[str], T]:
    def _c(x):
        try:
            return constructor(x)
        except PoagraphBuildException as p:
            raise argparse.ArgumentError(f"Incorrect argument {x}") from p
    return _c


def get_parser() -> argparse.ArgumentParser:
    """Create ArgumentParser for pang module."""

    p = argparse.ArgumentParser(prog='pang',
                                description='This software builds poagraph and generates consensuses.',
                                epilog='For more information check github.com/meoke/pang')
    p.add_argument('--multialignment', '-m',
                   type=_maf_file,
                   required=True,
                   help='Path to the multialignment file. ' + inspect.getdoc(Maf))
    p.add_argument('--datatype',
                   type=_data_type,
                   default=DataType.Nucleotides,
                   help='\'n\' for nucleotides, \'p\' for proteins. ' + inspect.getdoc(DataType))
    p.add_argument('--metadata',
                   type=_metadata_file,
                   help='Path to the csv file. ' + inspect.getdoc(MetadataCSV))
    p.add_argument('-raw_maf',
                   action='store_true',
                   default=False,
                   help='Poagraph building from maf file parameter.'
                        'Set if the maf content must not be transformed to DAG before building poagraph. '
                        'Poagraph that was build in this way provides consensuses tree but the consensuses do not '
                        'reflect the real life sequences.')
    p.add_argument('-fasta_complementation',
                   type=_fasta_provider_option,
                   default=FastaProviderOption.NCBI,
                   help='\'ncbi\' for NCBI, \'file\' for file. ' + inspect.getdoc(FastaProviderOption))
    p.add_argument('-email',
                   type=cli_arg(EmailAddress),
                   help=inspect.getdoc(EmailAddress))

#     p.add_argument('-cache',
#                    action='store_true',
#                    help='Used if Fasta Complementation Option is \"NCBI\" '
#                         'Stores sequences downloaded from NCBI on local disc.'
#                         'They are reused between uses of this program.')
#     p.add_argument('-missing_n',
#                    type=str,
#                    help='If fasta_complementation is NO, a custom symbol for missing nucleotides can be specified.'
#                         'Make sure it is included in BLOSUM matrix you use.')
#     p.add_argument('--fasta_source_file', '-f',
#                    type=_file_arg,
#                    help='ZIP archive with fasta files used to complement missing parts of sequences in maf file.')
#
#
#
#     p.add_argument('--blosum',
#                    type=_file_arg,
#                    help='Path to the BLOSUM matrix used in consensus generation algorithm.'
#                         'If fasta_complementation option is NO and a custom symbol is provided, '
#                         'the matrix specified here must include this symbol.'
#                         'If fasta_complementation option is NO and a custom symbol is not provided, '
#                         'the matrix specified here must include symbol \'?\' '
#                         'as this is the default symbol for missing nucleotide.'
#                    )
#     p.add_argument('--output', '-o',
#                    type=_dir_arg,
#                    default=create_default_output_dir(Path(getcwd())),
#                    help='Output directory path.')
#     p.add_argument('-fasta',
#                    action='store_true',
#                    help='Set if fasta files for consensuses must be produced.')
#     p.add_argument('-output_po',
#                    action='store_true',
#                    default=False,
#                    help='Set if po file with entire pangraph (without any consensuses) must be produced.')
#     p.add_argument('-consensus',
#                    type=_consensus_algorithm_option,
#                    default=ConsensusAlgorithm.NO,
#                    help='Set if consensus must be generated. Values to choose: \'simple\' or \'tree\'.')
#     p.add_argument('-hbmin',
#                    type=_float_0_1,
#                    default=0.6,
#                    help='Simple POA algorithm parameter. '
#                         'The minimum value of sequence compatibility to generated consensus.')
#     p.add_argument('-r',
#                    nargs=2,
#                    type=_float_0_1,
#                    action=_RangeArgAction,
#                    default=[0, 1],
#                    help='Tree POA algorithm parameter.'
#                         'Specify what part of sorted capabilities should be searched for node cutoff. E.g. [0.2,0.8]')
#     p.add_argument('-multiplier',
#                    type=float,
#                    default=1,
#                    help='Tree POA algorithm parameter.'
#                         'Cutoff value for node parameter. The greater it is, the more granular the tree is.')
#     p.add_argument('-stop',
#                    type=_float_0_1,
#                    default=0.99,
#                    help='Tree POA algorithm parameter.'
#                         'Value of node compatibility above which the node is no more split.')
#     p.add_argument('-re_consensus',
#                    action='store_true',
#                    default=False,
#                    help='Tree POA algorithm parameter.'
#                         'Set if after producing children nodes, sequences should be moved to'
#                         ' siblings nodes if compatibility to its consensus is higher.')
# c
#     p.add_argument('-p',
#                    type=float,
#                    default=1,
#                    help='Tree consensus algorithm parameter.'
#                         'When deciding about consensus node split, the compatibilities are raised to the power o p.'
#                         'It enables to change the linear meaning of compatibility values.'
#                         'For p from range [0,1] it decreases distances between small compatibilities and '
#                         'increases distances between the bigger ones.'
#                         'For p > 1 it increases distances between small compatibilities and '
#                         'decreases distances between the bigger ones.'
#                    )
#     p.add_argument('-max',
#                    default=MaxCutoffOption.MAX2,
#                    type=_max_cutoff_option,
#                    help='Specify which strategy - MAX1 or MAX2 use '
#                         'for finding max cutoff (see details in README.md)')
#     p.add_argument('-node',
#                    default=NodeCutoffOption.NODE3,
#                    type=_node_cutoff_option,
#                    help='Specify which strategy - NODE1 (1), NODE2 (2), NODE3 (3) or NODE4 (4) use '
#                         'for finding max cutoff (see details in README.md)')
#     p.add_argument('-v', '--verbose',
#                    action='store_true',
#                    default=False,
#                    help='Set if detailed log files must be produced.')
#     p.add_argument('-q', '--quiet',
#                    action='store_true',
#                    default=False,
#                    help='Set to turn off console logging .')
#     p.add_argument('-output_with_nodes',
#                    action='store_true',
#                    default=False,
#                    help='Set if output json should include nodes (it significantly increases file size).')
    return p
