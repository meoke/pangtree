from .SequenceMetadata import SequenceMetadata
from typing import Dict

GenomeMetadataDict = Dict[int, SequenceMetadata]


class MultialignmentMetadata:
    def __init__(self, title: str, source: str, genomes_metadata: GenomeMetadataDict):
        self.title = title
        self.source = source
        self.genomes_metadata = genomes_metadata
