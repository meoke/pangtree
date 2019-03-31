from typing import NewType

SequenceID = NewType('SequenceID', str)

ColumnID = NewType('ColumnID', int)

NodeID = NewType('NodeID', int)

BlockID = NewType('BlockID', int)

Sequence = NewType('Sequence', str)

Base = NewType('Base', bytes)


def make_base(x: str) -> Base:
    if len(x) == 0:
        raise Exception("Empty string for base.")
    if len(x) > 1:
        raise Exception("Base string mus thave length 1.")
    return Base(str.encode(x))

