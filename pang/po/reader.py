from pathlib import Path

from graph.Pangraph import Pangraph
from metadata import MultialignmentMetadata


def read(path: Path, genomes_info: MultialignmentMetadata) -> Pangraph:
    print("czytam")
