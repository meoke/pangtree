from pathlib import Path

from consensus.input_types import Blosum, Hbmin


def get_simple_consensus_tree(blosum: Blosum, output_dir: Path, hbmin: Hbmin, verbose: bool):
    print(blosum.filepath)
    print(output_dir)
    print(hbmin.value)
    print(verbose)
