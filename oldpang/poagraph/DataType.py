from enum import Enum


class DataType(Enum):
    Nucleotides = 0
    Proteins = 1

    @staticmethod
    def get_parameter_description():
        return "TODO"
