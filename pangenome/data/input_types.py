from io import StringIO
from typing import Optional
from pathlib import Path


class InputError(Exception):
    pass


class Maf:
    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        self.filename = filename
        self._raise_exception_if_incorrect(filecontent)
        self.filecontent = filecontent

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if False:
            raise InputError("Incorrect maf because...")

    @staticmethod
    def get_parameter_description():
        return 'Accepted formats: .maf, .po.'


class MetadataCSV:
    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        self.filename = filename
        self._raise_exception_if_incorrect(filecontent)
        self.filecontent = filecontent

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if False:
            raise InputError("Incorrect maf because...")

    @staticmethod
    def get_descriprion():
        return  "See... examples\\Ebola\\ebola_metadata.csv"

    def get_all_sequences_ids(self):
        return ["seq1"]

    def get_sequence_metadata(self, seq_id):
        return {}


class Po:
    def __init__(self, filecontent: StringIO, filename: Optional[Path]):
        self.filename = filename
        self._raise_exception_if_incorrect(filecontent)
        self.filecontent = filecontent

    @staticmethod
    def _raise_exception_if_incorrect(filecontent):
        if False:
            raise InputError("Incorrect maf because...")


class MissingSymbol:
    def __init__(self, symbol: str) -> 'MissingSymbol':
        if len(symbol) != 1:
            raise InputError('Missing symbol must be a single character.')
