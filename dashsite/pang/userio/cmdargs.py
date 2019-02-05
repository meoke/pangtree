import argparse
from os import getcwd
from pathlib import Path
from typing import Union, Dict

from .pathtools import create_default_output_dir

ArgType = Union[str, float, str, Path]
ArgsList = Dict[str, ArgType]


def _file_arg(path: str) -> Path:
    """Check if path exists."""

    file_path = Path(path)
    if not file_path.is_file():
        raise argparse.ArgumentTypeError(f"File {path} does not exist or is not a file.")
    return file_path


def _dir_arg(path: str) -> Path:
    """Check if dir exists and creates it if not."""

    dir_path = Path(path)
    if not dir_path.exists():
        dir_path.mkdir()
    return dir_path


def _float_0_1(arg: str) -> float:
    """Check if convertable to float and in range [0,1]."""

    try:
        v = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{arg} was passed, a float excpected.")
    if v < 0 or v > 1:
        raise argparse.ArgumentTypeError(f"This argument must be in range [0,1].")
    return v


class _RangeArgAction(argparse.Action):
    """Command line argument \'range\' (\'-r\') validation"""

    def __call__(self, parser, namespace, values, option_string=None):
        if values[1] < values[0]:
            raise ValueError("First r argument must be smaller or equal than the second r argument")
        setattr(namespace, self.dest, values)


def _get_parser():
    """Create ArgumentParser for pang module."""

    p = argparse.ArgumentParser(prog='pang',
                                description='Consensus generation and visulization of Pangenome',
                                epilog='For more information check github.com/meoke/pang')
    p.add_argument('--multialignment', '-m',
                   type=_file_arg,
                   required=True,
                   help='Path to the mulitalignment file. Accepted formats: .maf, .po.')
    p.add_argument('--data', '-d',
                   type=_file_arg,
                   required=True,
                   help='Path to the json file with genomes specification. See... examples\Ebola\ebola_metadata.json')
    p.add_argument('--output', '-o',
                   type=_dir_arg,
                   default=create_default_output_dir(Path(getcwd())),
                   help='Output directory path.')
    p.add_argument('-fasta',
                   action='store_true',
                   help='Set if fasta files must be produced.')
    p.add_argument('-vis',
                   action='store_true',
                   help='Set if visualization must be produced.')
    p.add_argument('-consensus',
                   choices=['simple', 'tree'],
                   help='Set if consensus must be generated. Algorithms to choose: \'simple\' or \'tree\'.')
    p.add_argument('-hbmin',
                   type=_float_0_1,
                   default=0.6,
                   help='Simple POA algorithm parameter. '
                        'The minimum value of sequence compatibility to generated consensus')
    p.add_argument('-r',
                   nargs=2,
                   type=_float_0_1,
                   action=_RangeArgAction,
                   default=[0, 1],
                   help='Tree POA algorithm parameter.'
                        'Specify what part of sorted capabilities should be searched for node cutoff. E.g. [0.2,0.8]')
    p.add_argument('-multiplier',
                   type=float,
                   default=1,
                   help='Tree POA algorithm parameter.'
                        'Cutoff value for node parameter. The greater it is, the more granular the tree is.')
    p.add_argument('-stop',
                   type=_float_0_1,
                   default=0.99,
                   help='Tree POA algorithm parameter.'
                        'Value of node compatibility above which the node is no more split.')
    p.add_argument('-re_consensus',
                   action='store_true',
                   default=True,
                   help='Tree POA algorithm parameter.'
                        'Set if after producing children nodes, sequences should be moved to'
                        ' siblings nodes if compatibility to its consensus is higher.')
    p.add_argument('-anti_granular',
                   action='store_true',
                   default=True,
                   help='Tree POA algorithm parameter.'
                        'Set if consensuses tree should be processed in a way to avoid fragmentation.')
    p.add_argument('-fasta_complementation',
                   nargs='*',
                   help='Maf file usually contains not full sequences but only parts of them, aligned to each other. '
                        'To build an exact pangraph the full sequences must be retrieved from: '
                        'ncbi or local file system. '
                        'Don\'t use this argument if you want the pangraph to be build without full sequences.'
                        'Use it without additional argument if you want to download the lacking fragments from ncbi'
                        '(then make sure that sequence identifiers used in maf match the ncbi accession identifiers) '
                        'Use it with additional argument (path to the directory with fasta files) if you want '
                        'to use fasta from local file system (make sure sequence identifiers in maf match file names')
    return p


def get_validated_args():
    """Parse and validate command line arguments"""

    parser = _get_parser()
    try:
        return parser.parse_args()
    except Exception as e:
        raise parser.error(e)


def get_fasta_complementation_option(fasta_complementation):
    if fasta_complementation is None:
        return None
    elif fasta_complementation == []:
        return 'ncbi'
    elif fasta_complementation[0]:
        dir_path = Path(fasta_complementation[0])
        if not dir_path.exists():
            raise Exception("Fasta directory does not exist.")
        return dir_path
