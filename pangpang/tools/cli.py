import argparse
import inspect
import os
from io import StringIO
from pathlib import Path
from typing import TypeVar, Callable, Optional, Union, List

from consensus.cutoffs import FindMaxCutoff, MAX2, MAX1, NODE3, FindCutoff, FindNodeCutoff, NODE1, NODE2, NODE4
from consensus.input_types import Blosum, Hbmin, Range
from datamodel.DataType import DataType
from datamodel.builders import PoagraphBuildException
from datamodel.fasta_providers.FastaProvider import FastaProvider, UseCache
from datamodel.fasta_providers.FromNCBI import EmailAddress
from datamodel.input_types import Maf, MetadataCSV, Po, MissingSymbol

from datamodel.fasta_providers import FastaProvider
from datamodel.fasta_providers.ConstSymbolProvider import ConstSymbolProvider
from datamodel.fasta_providers.FromNCBI import FromNCBI
from datamodel.fasta_providers.FromFile import FromFile

from consensus import input_types as consensus_input_types

from tools import pathtools


class InvalidPath(Exception):
    pass


def _get_file_extension(arg: str) -> str:
    """Returns file extension if it is included in arg, else throws InvalidPath."""

    file_path_suffix = Path(arg).suffix
    try:
        return file_path_suffix.split('.')[1]
    except Exception:
        raise InvalidPath(f"Cannot find file extension in {arg}.")



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


# def _fasta_provider_option(arg_fasta_provider_option: str) -> FastaProviderOption:
#     """Converts command line argument to FastaProviderOption"""
#
#     try:
#         return FastaProviderOption[arg_fasta_provider_option.upper()]
#     except KeyError:
#         raise argparse.ArgumentError("Incorrect FASTA_PROVIDER argument.")


def _path_if_valid(path: str) -> Path:
    """Check if path exists."""

    file_path = Path(path)
    if not pathtools.file_exists(file_path):
        raise InvalidPath(f"{file_path}")
    return file_path


def _cli_dir_arg(path: str) -> Path:
    """Check if dir exists and creates it if not."""

    dir_path = Path(path)
    if not pathtools.dir_exists(dir_path):
        pathtools.create_dir(dir_path)
    return dir_path


T = TypeVar('T')


def _cli_file_arg(arg: str, constructor: Callable[[StringIO, Optional[Path]], T]) -> T:
    try:
        filepath = _path_if_valid(arg)
    except InvalidPath:
        raise argparse.ArgumentTypeError(f"File {arg} does not exist or is not a file.")
    with open(arg) as infile:
        file_content = StringIO(infile.read())
        try:
            return constructor(file_content, filepath)
        except PoagraphBuildException as p:
            raise argparse.ArgumentError("Incorrect file content") from p


def _mulitalignment_file(x: str) -> Union[Maf, Po]:
    file_type = _get_file_extension(x)
    if file_type == 'maf':
        return _cli_file_arg(x, Maf)
    elif file_type == 'po':
        return _cli_file_arg(x, Po)
    else:
        raise InvalidPath("Only multialignment files with .maf or .po can be processed.")


def _metadata_file(x: str) -> MetadataCSV:
    return _cli_file_arg(x, MetadataCSV)


def _blosum_file(x: str) -> Blosum:
    return _cli_file_arg(x, Blosum)


def _range_arg(x: List[str]) -> Range:
    return Range(x)


def _cli_arg(constructor: Callable[[str], T]) -> Callable[[str], T]:
    def _c(x):
        try:
            return constructor(x)
        except PoagraphBuildException as p:
            raise argparse.ArgumentError(f"Incorrect argument {x}") from p
    return _c

class _RangeArgAction(argparse.Action):
    """Command line argument \'range\' (\'-r\') validation"""

    def __call__(self, parser, namespace, values, option_string=None):
        if values[1] < values[0]:
            raise ValueError("First r argument must be smaller or equal than the second r argument")
        setattr(namespace, self.dest, values)


def get_parser() -> argparse.ArgumentParser:
    """Create ArgumentParser for pang module."""

    p = argparse.ArgumentParser(prog='pang',
                                description='This software builds poagraph and generates consensuses.',
                                epilog='For more information check github.com/meoke/pang')
    p.add_argument('-output_dir', '-o',
                   type=_cli_dir_arg,
                   default=get_default_output_dir(),
                   help='Output directory path.')
    p.add_argument('--multialignment', '-m',
                   type=_mulitalignment_file,
                   required=True,
                   help='Path to the multialignment file.')
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
    p.add_argument('-fasta_provider',
                   # type=_fasta_provider_option,
                   choices=['ncbi', 'file'],
                   help='\'ncbi\' for NCBI, \'file\' for file. MISSING_SYMBOL will be used if not set. ')
    p.add_argument('-missing_symbol',
                   metavar='MISSING_SYMBOL',
                   type=_cli_arg(MissingSymbol),
                   default=MissingSymbol(),
                   help=inspect.getdoc(MissingSymbol))
    p.add_argument('-email',
                   type=_cli_arg(EmailAddress),
                   help=inspect.getdoc(EmailAddress))
    p.add_argument('-cache',
                   action='store_true',
                   help='Set if fastas downloaded from ncbi should be cached locally in .fastacache folder. '
                        + inspect.getdoc(UseCache))
    p.add_argument('--fasta_file', '-f',
                   type=_path_if_valid,
                   help='ZIP archive with fasta files or fasta file used as missing nucleotides/proteins source.')
    p.add_argument('-consensus',
                   choices=['poa', 'tree'],
                   help='\'poa\' for direct result of poa software, \'tree\' for Consensuses Tree algorith.')
    p.add_argument('-blosum',
                   type=_blosum_file,
                   help='Path to the blosum file. ' + inspect.getdoc(Blosum))
    p.add_argument('-hbmin',
                   type=_cli_arg(Hbmin),
                   default=consensus_input_types.Hbmin(),
                   help='Simple POA algorithm parameter. '
                        'Hbmin value. ' + inspect.getdoc(Hbmin))
    p.add_argument('-max',
                   default='max2',
                   choices=['max1', 'max2'],
                   help='Tree POA algorithm parameter. ' +
                        'Specify which strategy - MAX1 or MAX2 use for finding max cutoff.')
    p.add_argument('-node',
                   default='node3',
                   choices=['node1', 'node2', 'node3', 'node4'],
                   help='Tree POA algorithm parameter. ' +
                        'Specify which strategy - NODE1, NODE2, NODE3 or NODE4 use for finding max cutoff.')
    p.add_argument('-r',
                   nargs=2,
                   action='append',
                   help='Tree POA algorithm, MAX1 strategy parameter. ' + inspect.getdoc(Range))
    p.add_argument('-multiplier',
                   type=_cli_arg(consensus_input_types.Multiplier),
                   default=consensus_input_types.Multiplier(),
                   help='Tree POA algorithm parameter.' + inspect.getdoc(consensus_input_types.Multiplier))
    p.add_argument('-stop',
                   type=_cli_arg(consensus_input_types.Stop),
                   default=consensus_input_types.Stop(),
                   help='Tree POA algorithm parameter.' + inspect.getdoc(consensus_input_types.Stop))
    p.add_argument('-p',
                   type=_cli_arg(consensus_input_types.P),
                   default=consensus_input_types.P(),
                   help='Tree consensus algorithm parameter.' + inspect.getdoc(consensus_input_types.P))
#     p.add_argument('-fasta',
#                    action='store_true',
#                    help='Set if fasta files for consensuses must be produced.')
    p.add_argument('-output_po',
                   action='store_true',
                   default=False,
                   help='Set if output must containt poagraph as .po file.')
#     p.add_argument('-consensus',
#                    type=_consensus_algorithm_option,
#                    default=ConsensusAlgorithm.NO,
#                    help='Set if consensus must be generated. Values to choose: \'simple\' or \'tree\'.')




#     p.add_argument('-node',
#                    default=NodeCutoffOption.NODE3,
#                    type=_node_cutoff_option,
#                    help='Specify which strategy - NODE1 (1), NODE2 (2), NODE3 (3) or NODE4 (4) use '
#                         'for finding max cutoff (see details in README.md)')
    p.add_argument('-v', '--verbose',
                   action='store_true',
                   default=False,
                   help='Set if detailed log files must be produced.')
#     p.add_argument('-q', '--quiet',
#                    action='store_true',
#                    default=False,
#                    help='Set to turn off console logging .')
#     p.add_argument('-output_with_nodes',
#                    action='store_true',
#                    default=False,
#                    help='Set if output json should include nodes (it significantly increases file size).')
    return p


def resolve_fasta_provider(args: argparse.Namespace) -> FastaProvider:
    if args.fasta_provider is None:
        return ConstSymbolProvider(args.missing_symbol)
    elif args.fasta_provider == 'ncbi':
        if args.email is None:
            raise Exception("Email address must be specified. It must be provided when fasta source is \'ncbi\'.")
        use_cache = args.cache if args.cache else False
        return FromNCBI(args.email, use_cache)
    elif args.fasta_provider == 'file':
        if args.fasta_file is None:
            raise Exception("Fasta file source must be specified. It must be provided when fasta source is \'local\'.")
        return FromFile(args.fasta_file)
    else:
        raise Exception("Not known fasta provider."
                        "Should be \'ncbi\' or \'file\' or None."
                        "Cannot build pangraph.")


def resolve_max_strategy(args: argparse.Namespace) -> FindMaxCutoff:
    if args.max is None or args.max == "max2":
        return MAX2()
    elif args.max == "max1":
        r = consensus_input_types.Range(args.r)
        return MAX1(r)
    else:
        raise Exception("Not known max cutoff strategy."
                        "Should be \'max1\' or \'max2\' or None."
                        "Cannot generate Consensus Tree.")


def resolve_node_strategy(args: argparse.Namespace) -> FindNodeCutoff:
    if args.node is None or args.node == "node3":
        return NODE3()
    elif args.node == "node1":
        return NODE1(args.multiplier)
    elif args.node == "node2":
        return NODE2(args.r)
    elif args.node == "node4":
        return NODE4()
    else:
        raise Exception("Not known node cutoff strategy."
                        "Should be \'node1\', \'node2\',\'node3\', \'node3\' or None."
                        "Cannot generate Consensus Tree.")


def get_default_output_dir():
    """Creates timestamped child dir under current working directory."""

    current_dir = pathtools.get_cwd()
    current_time = pathtools.get_current_time()
    output_dir_name = "_".join(["output/", current_time])
    output_dir_path = pathtools.get_child_path(current_dir, output_dir_name)
    pathtools.create_dir(output_dir_path)
    return output_dir_path


def get_default_blosum(missing_base_symbol: MissingSymbol):
    """Returns default blosum file: Blosum80.mat"""
    parent_dir = Path(os.path.dirname(os.path.abspath(__file__)) + '/')
    default_blosum_path = pathtools.get_child_path(parent_dir, "../../bin/blosum80.mat")
    blosum_content = pathtools.get_file_content_stringio(default_blosum_path)
    return Blosum(blosum_content, default_blosum_path, missing_base_symbol)
