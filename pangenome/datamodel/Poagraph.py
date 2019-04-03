from typing import Optional, List

from datamodel.DataType import DataType
from datamodel.Sequence import Sequences
from datamodel.Node import Node
from datamodel.builders import maf2poagraph, maf2dagmaf, dagmaf2poagraph
from datamodel.input_types import Po, Maf, MetadataCSV
from .fasta_providers.FastaProvider import FastaProvider


class Poagraph:
    def __init__(self,
                 nodes: List[Node],
                 sequences: Sequences,
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
    def build_from_dagmaf(cls,
                          maf: Maf,
                          fasta_provider: FastaProvider,
                          metadata: Optional[MetadataCSV],
                          datatype: Optional[DataType] = DataType.Nucleotides) -> 'Poagraph':
        dagmaf = maf2dagmaf.get_dagmaf(maf)
        nodes, sequences = dagmaf2poagraph.get_poagraph(dagmaf, fasta_provider, metadata)
        poagraph = Poagraph(nodes, sequences)
        if datatype is not None:
            poagraph.datatype = datatype
        return poagraph

    @classmethod
    def build_from_po(cls, po: Po) -> 'Poagraph':
        return cls("nodes po", "paths po")

    def __eq__(self, other: 'Poagraph') -> bool:
        return self.nodes == other.nodes and \
            self.sequences == other.sequences and \
            self.datatype == other.datatype


