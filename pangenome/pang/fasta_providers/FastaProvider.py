import abc


class FastaProvider(abc.ABC):
    @abc.abstractmethod
    def get_source(self, sequenceID: str, start: int = None, end: int = None):
        pass

    def get_raw_sequence_from_fasta(self, fasta_handle):
        _ = fasta_handle.readline()
        return fasta_handle.read().replace('\n', '')
