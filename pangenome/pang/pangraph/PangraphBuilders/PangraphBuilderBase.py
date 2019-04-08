import abc

from pangenome.pang.metadata.MultialignmentMetadata import MultialignmentMetadata


class PangraphBuilderBase(abc.ABC):
    def __init__(self, genomes_info: MultialignmentMetadata):
        self.sequences_ids = genomes_info.get_all_sequences_ids()

    @abc.abstractmethod
    def build(self, input_content, pangraph):
        pass
