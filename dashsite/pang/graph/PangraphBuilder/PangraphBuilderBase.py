import abc

from metadata.MultialignmentMetadata import MultialignmentMetadata


class PangraphBuilderBase(abc.ABC):
    def __init__(self, genomes_info: MultialignmentMetadata):
        self.sequences_names = genomes_info.get_all_mafnames()

    @abc.abstractmethod
    def build(self, input, pangraph):
        pass
