from typing import NewType, List

from data.builders import PoagraphBuildException
from .Node import NodeID

SequencePath = NewType('SequncePath', List[NodeID])


class SequenceID:
    value: str

    def __init__(self, sequence_id: str):
        splitted = sequence_id.split('.')
        if len(splitted) > 1:
            self.value = ".".join(splitted[1:])
        elif len(splitted) == 1:
            self.value = splitted[0]
        else:
            raise PoagraphBuildException("Sequence ID cannot be empty.")


class SequenceMetadata:
    def __init__(self):
        pass


class Sequence:
    def __init__(self, seqid: SequenceID, paths: List[SequencePath], seqmetadata: SequenceMetadata):
        self.seqid: SequenceID = seqid
        self.paths: List[SequencePath] = paths
        self.seqmetadata: SequenceMetadata = seqmetadata
