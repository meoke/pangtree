from .SequenceID import SequenceID
from .SequenceMetadata import SequenceMetadata
from typing import Dict


class MultialignmentMetadata:
    def __init__(self, title: str, version: str, genomes_metadata: Dict[SequenceID, SequenceMetadata]):
        self.title = title
        self.version = version
        self.genomes_metadata = genomes_metadata

    # def get_id(self, sequence_name: str) -> int:
    #     return [seq_id for seq_id, data in self.genomes_metadata.items() if data.mafname == sequence_name][0]
    #
    # def get_group(self, sequence_name: str) -> str:
    #     return [data.group for seq_id, data in self.genomes_metadata.items() if data.mafname == sequence_name][0]
    #
    # def get_title(self, sequence_name: str) -> str:
    #     return [data.title for seq_id, data in self.genomes_metadata.items() if data.mafname == sequence_name][0]

    def get_all_mafnames(self):
        return [data.mafname for seq_id, data in self.genomes_metadata.items()]

    def feed_with_dagmaf_data(self, mafcontent):
        #todo
        pass
