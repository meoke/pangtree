from pangtreebuild.datamodel.Node import Base
from pangtreebuild.datamodel.Sequence import SequenceID
from pangtreebuild.datamodel.fasta_providers.FastaProvider import FastaProvider
from pangtreebuild.datamodel.input_types import MissingSymbol


class ConstSymbolProvider(FastaProvider):
    def __init__(self, missing_symbol: MissingSymbol):
        self.missing_symbol: Base = Base(missing_symbol.value)

    def get_base(self, sequence_id: SequenceID, i: int) -> Base:
        return self.missing_symbol
