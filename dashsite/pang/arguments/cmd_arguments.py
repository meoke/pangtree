import argparse
from os import getcwd
from pathlib import Path
from typing import Union, Dict

from tools import pathtools
from arguments.PangenomeParameters import FastaComplementationOption, PangenomeParameters, ConsensusAlgorithm, MaxCutoffOption, \
    NodeCutoffOption
from tools.pathtools import create_default_output_dir

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


def _fasta_complementation_option(arg_fasta_complementation: str) -> FastaComplementationOption:
    """Converts command line argument to FastaComplementationOption"""

    try:
        return FastaComplementationOption[arg_fasta_complementation.upper()]
    except KeyError:
        raise argparse.ArgumentError("Incorrect FASTA_COMPLEMENTATION argument.")


def _max_cutoff_option(max_cutoff_option: str) -> MaxCutoffOption:
    """Converts command line argument to MaxCutoffOption"""

    try:
        return MaxCutoffOption[max_cutoff_option.upper()]
    except KeyError:
        raise argparse.ArgumentError()


def _consensus_algorithm_option(arg_consensus_option: str) -> ConsensusAlgorithm:
    """Converts command line argument to ConsensusOption"""

    try:
        return ConsensusAlgorithm[arg_consensus_option.upper()]
    except KeyError:
        raise argparse.ArgumentError()


def _node_cutoff_option(node_cutoff_option: str) -> NodeCutoffOption:
    """Converts command line argument to NodeCutoffOption"""

    try:
        return NodeCutoffOption[node_cutoff_option.upper()]
    except KeyError:
        raise argparse.ArgumentError()


class _RangeArgAction(argparse.Action):
    """Command line argument \'range\' (\'-r\') validation"""

    def __call__(self, parser, namespace, values, option_string=None):
        if values[1] < values[0]:
            raise ValueError("First r argument must be smaller or equal than the second r argument")
        setattr(namespace, self.dest, values)


def _get_parser() -> argparse.ArgumentParser:
    """Create ArgumentParser for pang module."""

    p = argparse.ArgumentParser(prog='pang',
                                description='Build pangraph and generate consensuses',
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
                   help='Set if fasta files for consensuses must be produced.')
    p.add_argument('-consensus',
                   type=_consensus_algorithm_option,
                   default=ConsensusAlgorithm.NO,
                   help='Set if consensus must be generated. Values to choose: \'simple\' or \'tree\'.')
    p.add_argument('-hbmin',
                   type=_float_0_1,
                   default=0.6,
                   help='Simple POA algorithm parameter. '
                        'The minimum value of sequence compatibility to generated consensus.')
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
    p.add_argument('-not_dag',
                   action='store_true',
                   default=False,
                   help='Pangraph building from maf file parameter.'
                        'Set if the maf content must not be transformed to DAG when building pangraph. '
                        'Pangraph that was build in this way provides consensuses tree the consensuses do not '
                        'reflect the real life sequences.')
    p.add_argument('-fasta_complementation',
                   default=FastaComplementationOption.NCBI,
                   type=_fasta_complementation_option,
                   help='Pangraph building from maf file parameter. Ignored when -not_dag parameter is set.'
                        'Maf file usually contains not full sequences but only parts of them, aligned to each other. '
                        'To build an exact pangraph the full sequences must be retrieved from: '
                        'ncbi or local file system. '
                        'Don\'t use this argument if you want the pangraph to be build without full sequences.'
                        'Pass "ncbi" if you want to download the lacking fragments from ncbi'
                        'Pass "local" if you want to use fasta from local file system.')
    p.add_argument('-fasta_dir',
                   type=_dir_arg,
                   help='Local directory with fasta files used to complement missing parts of sequences in maf file.')
    p.add_argument('-p',
                   type=float,
                   default=1,
                   help='Tree consensus algorithm parameter.'
                        'When finding compatibilities cutoff, their values are raised to the power o p.')
    p.add_argument('-max',
                   default=MaxCutoffOption.MAX2,
                   type=_max_cutoff_option,
                   help='Specify which strategy - MAX1 or MAX2 use '
                        'for finding max cutoff (see details in README.md)')
    p.add_argument('-node',
                   default=NodeCutoffOption.NODE3,
                   type=_node_cutoff_option,
                   help='Specify which strategy - NODE1 (1), NODE2 (2), NODE3 (3) or NODE4 (4) use '
                        'for finding max cutoff (see details in README.md)')
    return p


def create_pangenome_parameters() -> PangenomeParameters:
    """Parse and validate command line arguments"""

    parser = _get_parser()
    args = parser.parse_args()
    return PangenomeParameters(
            multialignment_file_content=pathtools.get_file_content_as_stringio(args.multialignment),
            multialignment_file_path=args.multialignment,
            metadata_file_content=pathtools.get_file_content(args.data),
            metadata_file_path=args.data,
            output_path=args.output,
            generate_fasta=args.fasta,
            consensus_type=args.consensus,
            hbmin=args.hbmin,
            range=args.r,
            multiplier=args.multiplier,
            stop=args.stop,
            re_consensus=args.re_consensus,
            not_dag=args.not_dag,
            fasta_complementation_option=args.fasta_complementation,
            local_fasta_dirpath=args.fasta_dir,
            max_cutoff_option=args.max,
            node_cutoff_option=args.node
        )
