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
                   help='Path to the file with mulitalignment. Accepted formats: .maf, .po.')
    p.add_argument('--data', '-d',
                   type=_file_arg,
                   required=True,
                   help='Path to the json file with genomes specification. Format details: TODO.')
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
                   help='Use if consensus must be generated. Algorithms to choose: \'simple\' or \'tree\'.')
    p.add_argument('-hbmin',
                   type=_float_0_1,
                   default=0.5,
                   help='POA algorithm parameter. TODO')
    p.add_argument('-mincomp',
                   type=_float_0_1,
                   default=0.5,
                   help='Tree POA algorithm parameter. TODO')
    p.add_argument('-r',
                   nargs=2,
                   type=_float_0_1,
                   action=_RangeArgAction,
                   default=[0.2, 0.5],
                   help='Tree POA algorithm parameter. TODO')
    p.add_argument('-multiplier',
                   type=float,
                   default=1,
                   help='Tree POA algorithm parameter. TODO')
    p.add_argument('-stop',
                   type=_float_0_1,
                   default=0.99,
                   help='Tree POA algorithm parameter. TODO')
    p.add_argument('-re_consensus',
                   action='store_true',
                   default=True,
                   help='Tree POA algorithm parameter. TODO')
    return p


def get_validated_args():
    """Parse and validate command line arguments"""

    parser = _get_parser()
    try:
        return parser.parse_args()
    except Exception as e:
        raise parser.error(e)
