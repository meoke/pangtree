from typing import NewType

SequenceID = NewType('SequenceID', str)

ColumnID = NewType('ColumnID', int)

NodeID = NewType('NodeID', int)

BlockID = NewType('BlockID', int)

Sequence = NewType('Sequence', str)

Nucleobase = NewType('Nucleobase', bytes)


def make_nucleobase(x: str) -> Nucleobase:
    if len(x) == 0:
        raise Exception("Empty string for nucleobase.")
    if len(x) > 1:
        raise Exception("Nucleobase string mus thave length 1.")
    return Nucleobase(str.encode(x))

