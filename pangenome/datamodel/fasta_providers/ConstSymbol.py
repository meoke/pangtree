from .FastaProvider import FastaProvider
from ..input_types import MissingSymbol


class ConstSymbol(FastaProvider):
    def __init__(self, missing_symbol: MissingSymbol):
        self.missing_symbol = missing_symbol

    def get_base(self, sequence_id: str, i: int):
        return self.missing_symbol
