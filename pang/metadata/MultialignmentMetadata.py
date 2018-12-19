from .SequenceMetadata import SequenceMetadata
from typing import Dict


class MultialignmentMetadata:
    def __init__(self, title: str, version: str, genomes_metadata: Dict[int, SequenceMetadata]):
        self.title = title
        self.version = version
        self.genomes_metadata = genomes_metadata

    def get_id(self, sequence_name: str) -> int:
        return [seq_id for seq_id, data in self.genomes_metadata.items() if data.mafname == sequence_name][0]

    def get_group(self, sequence_name: str) -> str:
        a = [data.group for seq_id, data in self.genomes_metadata.items() if data.mafname == sequence_name][0]
        return a
