import argparse
import inspect
import os
from io import StringIO
from pathlib import Path
from typing import TypeVar, Callable, Optional, Union

from pangtreebuild.serialization.json import TaskParameters
from pangtreebuild.affinity_tree import parameters as at_params
from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.pangenome.parameters import missings
from pangtreebuild.tools import pathtools


class InvalidPath(Exception):
    pass


def _get_file_extension(arg: str) -> str:
    """Returns file extension if present in arg else throws InvalidPath."""

    file_path_suffix = Path(arg).suffix
    try:
        return file_path_suffix.split('.')[1]
    except Exception:
        raise InvalidPath(f"Cannot find file extension in {arg}.")


def _data_type(data_type: str) -> graph.DataType:
    """Converts command line argument to DataType"""

    if data_type == "p":
        data_type_full_name = "Proteins"
    elif data_type == "n":
        data_type_full_name = "Nucleotides"
    else:
        raise argparse.ArgumentError("""Unknown data type. \'p\' for proteins
                                        or \'n\' for nucleotides available.""")
    try:
        dt = graph.DataType[data_type_full_name]
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


def _cli_file_arg(arg: str,
                  constructor: Callable[[StringIO, Optional[Path]], T]) -> \
        T:
    try:
        filepath = _path_if_valid(arg)
    except InvalidPath:
        raise argparse.ArgumentTypeError(f"""File {arg} does not exist
                                             or is not a file.""")
    with open(arg) as infile:
        file_content = StringIO(infile.read())
        try:
            return constructor(file_content, filepath)
        except Exception as p:
            raise argparse.ArgumentError("Incorrect file content") from p


def _mulitalignment_file(x: str) -> Union[msa.Maf, msa.Po]:
    file_type = _get_file_extension(x)
    if file_type == 'maf':
        return _cli_file_arg(x, msa.Maf)
    elif file_type == 'po':
        return _cli_file_arg(x, msa.Po)
    else:
        raise InvalidPath("""Only msa files with
                            .maf or .po can be processed.""")


def _metadata_file(x: str) -> msa.MetadataCSV:
    return _cli_file_arg(x, msa.MetadataCSV)


def _blosum_file(x: str) -> at_params.Blosum:
    return _cli_file_arg(x, at_params.Blosum)


def _cli_arg(constructor: Callable[[str], T]) -> Callable[[str], T]:
    def _c(x):
        try:
            return constructor(x)
        except Exception as p:
            raise argparse.ArgumentError(f"Incorrect argument {x}") from p
    return _c


def get_parser() -> argparse.ArgumentParser:
    """Create ArgumentParser for pang module."""

    p = argparse.ArgumentParser(prog='pangtreebuild',
                                description="""This software builds poagraph
                                                and generates affinitytree.""",
                                epilog="""For more information check
                                          github.com/meoke/pangtree""")
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
                   default=graph.DataType.Nucleotides,
                   help='\'n\' for nucleotides, \'p\' for proteins. ' +
                        inspect.getdoc(graph.DataType))
    p.add_argument('--metadata',
                   metavar='METADATA_PATH',
                   type=_metadata_file,
                   help='Path to the csv file with metadata. ' +
                        inspect.getdoc(msa.MetadataCSV))
    p.add_argument('--raw_maf',
                   action='store_true',
                   default=False,
                   help="""Poagraph building from maf file parameter. Set if
                           the maf content must not be transformed to DAG
                           before building graph. Poagraph that was build
                           in this way provides affinitytree tree but the
                           affinitytree do not reflect the real life
                           sequences.""")
    p.add_argument('--fasta_provider',
                   metavar="FASTA_PROVIDER",
                   choices=['ncbi', 'file'],
                   help="""'Maf file may not include full _sequences.
                            In such case an additional data source is needed.
                            Use \'ncbi\' for NCBI (activates CACHE option)
                            or \'file\' for file (then provide also
                            FASTA_PATH). MISSING_SYMBOL is used if this
                            argument is omitted.""")
    p.add_argument('--missing_symbol',
                   metavar='MISSING_SYMBOL',
                   type=_cli_arg(missings.MissingBase),
                   help=inspect.getdoc(missings.MissingBase))
    p.add_argument('--cache',
                   action='store_true',
                   help="""Set if fastas downloaded from NCBI should be cached
                           locally in .fastacache folder. Used if Fasta
                           Provider is NCBI. Sequences downloaded from NCBI
                           are stored and reused by this program.""")
    p.add_argument('--fasta_path',
                   metavar="FASTA_PATH",
                   type=_path_if_valid,
                   help="""ZIP archive with fasta files or fasta file used
                        as FASTA_PROVIDER.""")
    p.add_argument('--affinity',
                   choices=['poa', 'tree'],
                   help="""Generate affinity tree. Use \'poa\' for direct
                           result of poa software, \'tree\' for Affinity
                           Tree algorithm.""")
    p.add_argument('--blosum',
                   type=_blosum_file,
                   metavar='BLOSUM_PATH',
                   help='Path to the blosum file. ' +
                        inspect.getdoc(at_params.Blosum))
    p.add_argument('--hbmin',
                   type=_cli_arg(at_params.Hbmin),
                   default=at_params.Hbmin(),
                   help='Simple POA algorithm parameter. '
                        'Hbmin value. ' + inspect.getdoc(at_params.Hbmin))
    p.add_argument('--stop',
                   type=_cli_arg(at_params.Stop),
                   default=at_params.Stop(),
                   help='Tree POA algorithm parameter.' +
                        inspect.getdoc(at_params.Stop))
    p.add_argument('--p',
                   type=_cli_arg(at_params.P),
                   default=at_params.P(),
                   help='Tree consensus algorithm parameter.' +
                        inspect.getdoc(at_params.P))
    p.add_argument('--output_fasta',
                   action='store_true',
                   help="""Set if fasta files for _sequences and
                            affinitytree must be produced.""")
    p.add_argument('--output_po',
                   action='store_true',
                   default=False,
                   help='Set if po file for poagraph must be produced.'),
    p.add_argument('--output_full',
                   action='store_true',
                   default=False,
                   help='Set if the result pangenome.json should contain '
                        'list of nodes ids for sequences and consensuses'),
    p.add_argument('-v', '--verbose',
                   action='store_true',
                   default=False,
                   help='Set if detailed log files must be produced.')
    p.add_argument('-q', '--quiet',
                   action='store_true',
                   default=False,
                   help='Set to turn off console logging.')
    return p


def resolve_fasta_provider(args: argparse.Namespace) -> \
        missings.FastaProvider:
    """Returns fasta provider based on parsed arguments."""

    if args.fasta_provider is None:
        if args.missing_symbol is None:
            return missings.ConstBaseProvider(missings.MissingBase())
        else:
            return missings.ConstBaseProvider(args.missing_symbol)
    elif args.fasta_provider == 'ncbi':
        use_cache = args.cache if args.cache else False
        return missings.FromNCBI(use_cache)
    elif args.fasta_provider == 'file':
        if args.fasta_path is None:
            raise Exception("""Fasta file source must be specified.
                               It must be provided when fasta source
                               is \'local\'.""")
        return missings.FromFile(args.fasta_path)
    else:
        raise Exception("""Not known fasta provider.
                           Should be \'ncbi\' or \'file\' or None.
                           Cannot build graph.""")


def get_default_output_dir():
    """Creates timestamped child dir under current working directory."""

    current_dir = pathtools.get_cwd()
    output_dir = pathtools.get_child_path(current_dir, "output")
    pathtools.create_dir(output_dir)
    current_time = pathtools.get_current_time()
    output_dir_name = "_".join(["output", current_time])
    output_dir_path = pathtools.get_child_path(output_dir,
                                               output_dir_name)
    pathtools.create_dir(output_dir_path)
    return output_dir_path


def get_default_blosum():
    """Returns default blosum file: Blosum80.mat"""

    pangtreebuild_dir = Path(__file__).parent.parent
    default_blosum_path = pathtools.get_child_path(pangtreebuild_dir, "affinity_tree/bin/blosum80.mat")
    blosum_content = pathtools.get_file_content_stringio(default_blosum_path)
    return at_params.Blosum(blosum_content, default_blosum_path)


def get_task_parameters(args: argparse.Namespace, running_time) -> \
        TaskParameters:
    """Returns TaskParameters object based on parsed arguments."""

    return TaskParameters(running_time=running_time,
                          multialignment_file_path=args.multialignment.filename,
                          multialignment_format=str(type(args.multialignment).__name__),
                          datatype=args.datatype.name,
                          metadata_file_path=args.metadata.filename if args.metadata else None,
                          blosum_file_path=args.blosum.filepath if args.blosum else None,
                          output_path=args.output_dir,
                          output_po=bool(args.output_po),
                          output_fasta=bool(args.output_fasta),
                          output_with_nodes=bool(args.output_full),
                          verbose=bool(args.verbose),
                          raw_maf=bool(args.raw_maf),
                          fasta_provider=args.fasta_provider if args.fasta_provider else 'ConstSymbol',
                          cache=bool(args.cache),
                          missing_base_symbol=args.missing_symbol.value if args.missing_symbol else missings.MissingBase().value,
                          fasta_source_file=str(args.fasta_path),
                          consensus_type=args.affinity,
                          hbmin=args.hbmin.value if args.hbmin else None,
                          stop=args.stop.value if args.stop else None,
                          p=args.p.value if args.p else None)
