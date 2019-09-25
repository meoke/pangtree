from typing import Optional, List, Tuple, Dict
from time import time
from datetime import datetime
from pangtreebuild.datamodel.fasta_providers.ConstSymbolProvider import ConstSymbolProvider
from pangtreebuild.affinity_tree.structure import Compatibility
from pangtreebuild.affinity_tree.parameters import P
from pangtreebuild.datamodel.DAGMaf import DAGMaf
from pangtreebuild.datamodel.DataType import DataType
from pangtreebuild.datamodel.Sequence import SeqPath, SequenceID, Sequence
from pangtreebuild.datamodel.Node import Node
from pangtreebuild.datamodel.builders import maf2poagraph, maf2dagmaf, dagmaf2poagraph, po2poagraph
from pangtreebuild.datamodel.input_types import Po, Maf, MetadataCSV, MissingSymbol
from pangtreebuild.datamodel.fasta_providers.FastaProvider import FastaProvider
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
        if metadata:
            Poagraph._complement_metadata_for_sequences_absent_in_metadata_provided(poagraph, metadata)
        poagraph.datatype = datatype
        return poagraph

    @classmethod
    def build_from_dagmaf(cls,
                          maf: Maf,
                          fasta_provider: Optional[FastaProvider] = ConstSymbolProvider(MissingSymbol()),
                          metadata: Optional[MetadataCSV] = None,
                          datatype: Optional[DataType] = DataType.Nucleotides) -> Tuple['Poagraph', DAGMaf]:
        dagmaf = maf2dagmaf.get_dagmaf(maf)
        nodes, sequences = dagmaf2poagraph.get_poagraph(dagmaf, fasta_provider, metadata)
        poagraph = Poagraph(nodes, sequences)
        if metadata:
            Poagraph._complement_metadata_for_sequences_absent_in_metadata_provided(poagraph, metadata)
        poagraph.datatype = datatype
        return poagraph, dagmaf

    @classmethod
    def build_from_po(cls,
                      po: Po,
                      metadata: Optional[MetadataCSV] = None,
                      datatype: Optional[DataType] = DataType.Nucleotides) -> 'Poagraph':
        nodes, sequences = po2poagraph.get_poagraph(po, metadata)
        poagraph = Poagraph(nodes, sequences)
        if metadata:
            Poagraph._complement_metadata_for_sequences_absent_in_metadata_provided(poagraph, metadata)
        poagraph.datatype = datatype
        return poagraph

    def __eq__(self, other: 'Poagraph') -> bool:
        return self.nodes == other.nodes and \
            self.sequences == other.sequences and \
            self.datatype == other.datatype

    def get_compatibilities(self,
                            sequences_ids: List[SequenceID],
                            consensus_path: SeqPath,
                            p: Optional[P] = P(1)) -> Dict[SequenceID, Compatibility]:
        compatibilities = dict()
        for seq_id in sequences_ids:
            try:
                sequence_paths = self.sequences[seq_id].paths
            except KeyError:
                raise Exception("No sequence with given ID in poagraph.")
            if len(sequence_paths) == 1:
                sequence_path = sequence_paths[0]
            else:
                sequence_path = [node_id for path in sequence_paths for node_id in path]
            compatibilities[seq_id] = Compatibility(len(set(sequence_path).intersection(set(consensus_path))) /
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
        # return {path_key: 100 for path_key in unweighted_sources_weights.keys()}
        return normalized_sources_weights_dict

    def get_sequence_nodes_count(self, seq_id: SequenceID) -> int:
        """Return the sum of lengths of all paths with the seid."""
        if seq_id not in self.sequences:
            raise Exception("No sequence with given ID in poagraph.")
        return sum([len(path) for path in self.sequences[seq_id].paths])

    def get_sequences_ids(self) -> List[SequenceID]:
        """Returns sequences of all genome sequences in Poagraph."""

        return [*self.sequences.keys()]

    @staticmethod
    def _complement_metadata_for_sequences_absent_in_metadata_provided(poagraph: 'Poagraph', metadata: MetadataCSV):
        headers = metadata.get_metadata_keys()
        for seq in poagraph.sequences.values():
            for h in headers:
                if h not in seq.seqmetadata:
                    seq.seqmetadata[h] = None
