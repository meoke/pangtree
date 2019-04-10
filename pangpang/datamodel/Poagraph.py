from typing import Optional, List, Tuple, Dict

from consensus.ConsensusTree import CompatibilityToPath
from consensus.input_types import P
from datamodel.DAGMaf import DAGMaf
from datamodel.DataType import DataType
from datamodel.Sequence import SequencePath, SequenceID, Sequence
from datamodel.Node import Node
from datamodel.builders import maf2poagraph, maf2dagmaf, dagmaf2poagraph
from datamodel.input_types import Po, Maf, MetadataCSV
from .fasta_providers.FastaProvider import FastaProvider
import numpy as np

class Poagraph:
    def __init__(self,
                 nodes: List[Node],
                 sequences: Dict[SequenceID, Sequence],
                 datatype: Optional[DataType] = DataType.Nucleotides):
        self.nodes: List[Node] = nodes
        self.sequences: Dict[SequenceID, Sequence] = sequences
        self.datatype: DataType = datatype

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
                          datatype: Optional[DataType] = DataType.Nucleotides) -> Tuple['Poagraph', DAGMaf]:
        dagmaf = maf2dagmaf.get_dagmaf(maf)
        nodes, sequences = dagmaf2poagraph.get_poagraph(dagmaf, fasta_provider, metadata)
        poagraph = Poagraph(nodes, sequences)
        if datatype is not None:
            poagraph.datatype = datatype
        return poagraph, dagmaf

    @classmethod
    def build_from_po(cls, po: Po) -> 'Poagraph':
        return cls("nodes po", "paths po")

    def __eq__(self, other: 'Poagraph') -> bool:
        return self.nodes == other.nodes and \
            self.sequences == other.sequences and \
            self.datatype == other.datatype

    def get_compatibilities(self,
                            sequences_ids: List[SequenceID],
                            consensus_path: SequencePath,
                            p: Optional[P] = P(1)) -> Dict[SequenceID, CompatibilityToPath]:
        compatibilities = dict()
        for seq_id in sequences_ids:
            try:
                sequence_paths = self.sequences[seq_id].paths
            except KeyError:
                raise Exception("No sequence with given ID in pangraph.")
            if len(sequence_paths) == 1:
                sequence_path = sequence_paths[0]
            else:
                sequence_path = [node_id for path in sequence_paths for node_id in path]
            compatibilities[seq_id] = CompatibilityToPath(len(set(sequence_path).intersection(set(consensus_path))) /
                                                          len(sequence_path), p)
        return compatibilities

    def get_sequences_weights(self, sequences_ids: List[SequenceID]):
        if not sequences_ids:
            return dict()

        a = np.zeros(len(self.nodes), dtype=np.int)
        unweighted_sources_weights = {}
        for seq_id in sequences_ids:
            for path in self.sequences[seq_id].paths:
                a[path] += 1

        for seq_id in sequences_ids:
            sequence_nodes_ids = [node_id for path in self.sequences[seq_id].paths for node_id in path]
            unweighted_sources_weights[seq_id] = np.mean(a[sequence_nodes_ids]) if any(a[sequence_nodes_ids]) else 0

        max_weight = max(unweighted_sources_weights.values())
        min_weight = min(unweighted_sources_weights.values())
        diff_weight = max_weight - min_weight
        if diff_weight == 0:
            normalized_sources_weights_dict = {path_key: 100 for path_key in unweighted_sources_weights.keys()}
        else:
            normalized_sources_weights_dict = {path: int((weight - min_weight)/diff_weight*100)
                                               for path, weight in unweighted_sources_weights.items()}
        return normalized_sources_weights_dict

    def get_sequence_nodes_count(self, seq_id: SequenceID) -> int:
        """Return the sum of lengths of all paths with the seid."""
        if seq_id not in self.sequences:
            raise Exception("No sequence with given ID in pangraph.")
        return sum([len(path) for path in self.sequences[seq_id].paths])


    def get_sequences_ids(self) -> List[SequenceID]:
        """Returns sequences of all genome sequences in Poagraph."""

        return [*self.sequences.keys()]

