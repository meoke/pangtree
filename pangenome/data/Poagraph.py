from typing import Optional, List, Dict

from data.DataType import DataType
from data.Node import Node
from data.Sequence import SequenceID, Sequence
from data.builders import maf2poagraph
from .input_types import Po, Maf, MetadataCSV
from .DAGMaf import DAGMaf
from .fasta_providers.FastaProvider import FastaProvider



class Poagraph:
    def __init__(self,
                 nodes: List[Node],
                 sequences: Dict[SequenceID, Sequence],
                 datatype: Optional[DataType] = DataType.Nucleotides):
        self.nodes = nodes
        self.sequences = sequences
        self.datatype = datatype

    @classmethod
    def build_from_maf(cls,
                       maf: Maf,
                       metadata: Optional[MetadataCSV] = None,
                       datatype: Optional[DataType] = DataType.Nucleotides) -> 'Poagraph':
        nodes, sequences = maf2poagraph.get_poagraph(maf, metadata)
        poagraph = Poagraph(nodes, sequences)
        if datatype is not None:
            poagraph.datatype = datatype
        return poagraph

    @classmethod
    def build_from_dagmaf(cls, dagmaf: DAGMaf, fasta_provider: FastaProvider, metadata: Optional[MetadataCSV]) -> 'Poagraph':
        return cls('nodes dagmaf', 'paths dagaf')

    @classmethod
    def build_from_po(cls, po: Po) -> 'Poagraph':
        return cls("nodes po", "paths po")


