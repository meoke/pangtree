from typing import NewType, List, Any, Dict

from .builders.PoagraphBuildException import PoagraphBuildException
from .Node import NodeID

SequencePath = NewType('SequencePath', List[NodeID])
SequenceMetadata = NewType('SequenceMetadata', Dict[str, Any])


class SequenceID:
    value: str

    def __init__(self, sequence_id: str, skip_part_before_dot=True):
        if not skip_part_before_dot:
            self.value = sequence_id
        else:
            splitted = sequence_id.split('.')
            if len(splitted) > 1:
                self.value = ".".join(splitted[1:])
            elif len(splitted) == 1:
                self.value = splitted[0]
            else:
                raise PoagraphBuildException("Sequence ID cannot be empty.")

    def __str__(self):
        return f"{self.value}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other: 'SequenceID'):
        return other and self.value.__eq__(other.value)

    def __hash__(self):
        return hash(self.value)


class Sequence:
    def __init__(self, seqid: SequenceID, paths: List[SequencePath], seqmetadata: SequenceMetadata):
        self.seqid: SequenceID = seqid
        self.paths: List[SequencePath] = paths
        self.seqmetadata: SequenceMetadata = seqmetadata

    def __str__(self):
        return f"seqid: {self.seqid}, paths: {[str(p) for p in self.paths]}, metadata: {str(self.seqmetadata)}"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other: 'Sequence'):
        return other and \
                self.seqid == other.seqid and \
                self.paths == other.paths and \
                self.seqmetadata == other.seqmetadata


Sequences = NewType("Sequences", Dict[SequenceID, Sequence])
