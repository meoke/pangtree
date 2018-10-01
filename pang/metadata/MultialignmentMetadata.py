from .SequenceMetadata import SequenceMetadata
from typing import Dict

GenomeMetadataDict = Dict[int, SequenceMetadata]


class MultialignmentMetadata:
    def __init__(self, title: str, version: str, genomes_metadata: GenomeMetadataDict):
        self.title = title
        self.version = version
        self.genomes_metadata = genomes_metadata
