import argparse
import inspect
import os
from io import StringIO
from pathlib import Path
from typing import TypeVar, Callable, Optional, Union, List

from pangtreebuild.output.PangenomeJSON import TaskParameters
from pangtreebuild.affinity_tree.parameters import Blosum, Hbmin
from pangtreebuild.datamodel.DataType import DataType
from pangtreebuild.datamodel.builders import PoagraphBuildException
from pangtreebuild.datamodel.fasta_providers.FastaProvider import UseCache
from pangtreebuild.datamodel.input_types import Maf, MetadataCSV, Po, MissingSymbol

from pangtreebuild.datamodel.fasta_providers import FastaProvider
from pangtreebuild.datamodel.fasta_providers.ConstSymbolProvider import ConstSymbolProvider
from pangtreebuild.datamodel.fasta_providers.FromNCBI import FromNCBI
from pangtreebuild.datamodel.fasta_providers.FromFile import FromFile

from pangtreebuild.affinity_tree import parameters as consensus_input_types

from pangtreebuild.tools import pathtools


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


def _cli_arg(constructor: Callable[[str], T]) -> Callable[[str], T]:
    def _c(x):
        try:
            return constructor(x)
        except PoagraphBuildException as p:
            raise argparse.ArgumentError(f"Incorrect argument {x}") from p
    return _c


def get_parser() -> argparse.ArgumentParser:
    """Create ArgumentParser for pang module."""

    p = argparse.ArgumentParser(prog='pangtreebuild',
                                description='This software builds poagraph and generates affinitytree.',
                                epilog='For more information check github.com/meoke/pang')
    p.add_argument('--output_dir',
                   type=_cli_dir_arg,
                   default=get_default_output_dir(),
                   help='Output directory path.')
    p.add_argument('--multialignment',
                   metavar='MULTIALIGNMENT_PATH',
                   type=_mulitalignment_file,
                   required=True,
                   help='Path to the multialignment file.')
    p.add_argument('--datatype',
                   type=_data_type,
                   default=DataType.Nucleotides,
                   help='\'n\' for nucleotides, \'p\' for proteins. ' + inspect.getdoc(DataType))
    p.add_argument('--metadata',
                   metavar='METADATA_PATH',
                   type=_metadata_file,
                   help='Path to the csv file with metadata. ' + inspect.getdoc(MetadataCSV))
    p.add_argument('--raw_maf',
                   action='store_true',
                   default=False,
                   help='Poagraph building from maf file parameter.'
                        'Set if the maf content must not be transformed to DAG before building poagraph. '
                        'Poagraph that was build in this way provides affinitytree tree but the affinitytree do not '
                        'reflect the real life sequences.')
    p.add_argument('--fasta_provider',
                   metavar="FASTA_PROVIDER",
                   choices=['ncbi', 'file'],
                   help='Maf file may not include full sequences. '
                        'In such case an additional data source is needed. '
                        'Use \'ncbi\' for NCBI (then CACHE option is available) or \'file\' for file (then provide also FASTA_PATH). MISSING_SYMBOL is used if this argument is omitted. ')
    p.add_argument('--missing_symbol',
                   metavar='MISSING_SYMBOL',
                   type=_cli_arg(MissingSymbol),
                   help=inspect.getdoc(MissingSymbol))
    p.add_argument('--cache',
                   action='store_true',
                   help='Set if fastas downloaded from NCBI should be cached locally in .fastacache folder. '
                        + inspect.getdoc(UseCache))
    p.add_argument('--fasta_path',
                   metavar="FASTA_PATH",
                   type=_path_if_valid,
                   help='ZIP archive with fasta files or fasta file used as FASTA_PROVIDER.')
    p.add_argument('--consensus',
                   choices=['poa', 'tree'],
                   help='Generate consensus tree. Use \'poa\' for direct result of poa software, \'tree\' for Consensus Tree algorithm.')
    p.add_argument('--blosum',
                   type=_blosum_file,
                   metavar='BLOSUM_PATH',
                   help='Path to the blosum file. ' + inspect.getdoc(Blosum))
    p.add_argument('--hbmin',
                   type=_cli_arg(Hbmin),
                   default=consensus_input_types.Hbmin(),
                   help='Simple POA algorithm parameter. '
                        'Hbmin value. ' + inspect.getdoc(Hbmin))
    p.add_argument('--stop',
                   type=_cli_arg(consensus_input_types.Stop),
                   default=consensus_input_types.Stop(),
                   help='Tree POA algorithm parameter.' + inspect.getdoc(consensus_input_types.Stop))
    p.add_argument('--p',
                   type=_cli_arg(consensus_input_types.P),
                   default=consensus_input_types.P(),
                   help='Tree consensus algorithm parameter.' + inspect.getdoc(consensus_input_types.P))
    p.add_argument('--output_fasta',
                   action='store_true',
                   help='Set if fasta files for sequences and affinitytree must be produced.')
    p.add_argument('--output_po',
                   action='store_true',
                   default=False,
                   help='Set if po file for poagraph must be produced.')
    p.add_argument('-v', '--verbose',
                   action='store_true',
                   default=False,
                   help='Set if detailed log files must be produced.')
    p.add_argument('-q', '--quiet',
                   action='store_true',
                   default=False,
                   help='Set to turn off console logging.')
    return p


def resolve_fasta_provider(args: argparse.Namespace) -> FastaProvider:
    if args.fasta_provider is None:
        if args.missing_symbol is None:
            return ConstSymbolProvider(MissingSymbol())
        else:
            return ConstSymbolProvider(args.missing_symbol)
    elif args.fasta_provider == 'ncbi':
        use_cache = args.cache if args.cache else False
        return FromNCBI(use_cache)
    elif args.fasta_provider == 'file':
        if args.fasta_path is None:
            raise Exception("Fasta file source must be specified. It must be provided when fasta source is \'local\'.")
        return FromFile(args.fasta_path)
    else:
        raise Exception("Not known fasta provider."
                        "Should be \'ncbi\' or \'file\' or None."
                        "Cannot build poagraph.")


def get_default_output_dir():
    """Creates timestamped child dir under current working directory."""

    current_dir = pathtools.get_cwd()
    output_dir = pathtools.get_child_path(current_dir, "output")
    pathtools.create_dir(output_dir)
    current_time = pathtools.get_current_time()
    output_dir_name = "_".join(["output", current_time])
    output_dir_path = pathtools.get_child_path(output_dir, output_dir_name)
    pathtools.create_dir(output_dir_path)
    return output_dir_path


def get_default_blosum():
    """Returns default blosum file: Blosum80.mat"""
    parent_dir = Path(os.path.dirname(os.path.abspath(__file__)) + '/')
    default_blosum_path = pathtools.get_child_path(parent_dir, "../../bin/blosum80.mat")
    blosum_content = pathtools.get_file_content_stringio(default_blosum_path)
    return Blosum(blosum_content, default_blosum_path)


def get_task_parameters(args: argparse.Namespace, running_time) -> TaskParameters:
    return TaskParameters(running_time=running_time,
                          multialignment_file_path=args.multialignment.filename,
                          multialignment_format=str(type(args.multialignment).__name__),
                          datatype=args.datatype.name,
                          metadata_file_path=args.metadata.filename if args.metadata else None,
                          blosum_file_path=args.blosum.filepath if args.blosum else None,
                          output_path=args.output_dir,
                          output_po=bool(args.output_po),
                          output_fasta=bool(args.output_fasta),
                          output_with_nodes=True,
                          verbose=bool(args.verbose),
                          raw_maf=bool(args.raw_maf),
                          fasta_provider=args.fasta_provider if args.fasta_provider else 'ConstSymbol',
                          cache=bool(args.cache),
                          missing_base_symbol=args.missing_symbol.value if args.missing_symbol else MissingSymbol().value,
                          fasta_source_file=args.fasta_path,
                          consensus_type=args.consensus,
                          hbmin=args.hbmin.value if args.hbmin else None,
                          stop=args.stop.value if args.stop else None,
                          p=args.p.value if args.p else None)