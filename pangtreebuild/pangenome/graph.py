"""Data model of poagraph"""

from enum import Enum
import numpy as np
from typing import Any, Dict, List, NewType, Optional, Union

from pangtreebuild.affinity_tree import parameters as at_params
from pangtreebuild.pangenome.parameters import msa


class DataType(Enum):
    """It describes whether poagraph is build from nucleotides or proteins."""

    Nucleotides = 0
    Proteins = 1


class Base(object):
    """Nucleotide or protein base.

    Args:
        base: base value as str

    Raises:
        ValueError: If base has length other than 1.

    Attributes:
        value (bytes): base value
    """

    def __init__(self, base: str):
        if len(base) == 0:
            raise ValueError("Poagraph node base cannot be empty.")
        if len(base) > 1:
            raise ValueError("Poagraph node base must have length = 1.")
        self.value: bytes = str.encode(base)

    def __eq__(self, other: 'Base'):
        return other and self.value.__eq__(other.value)

    def as_str(self) -> str:
        """Returns base value as string."""

        return self.value.decode("ASCII")


NodeID = NewType('NodeID', int)
ColumnID = NewType('ColumnID', int)
BlockID = NewType('BlockID', int)


class Node(object):
    """Poagraph node.

    Args:
        node_id: Node ID.
        base: Node base (single letter).
        aligned_to: ID of the poagraph node it is aligned to.
        column_id: Poagraph is visualised as succeeding columns of nodes,
                   so nodes have assigned a (not unique) column ID.
        block_id: If multialignment source was MAF this is ID of the block that the node is from. Otherwise == 1.

    Attributes:
        node_id: Node ID.
        base: Node base (single letter).
        aligned_to: ID of the poagraph node it is aligned to.
        column_id: Poagraph is visualised as succeeding columns of nodes,
                   so nodes have assigned a (not unique) column ID.
        block_id: If multialignment source was MAF this is ID of the block that the node is from. Otherwise == 1.
    """
    def __init__(self,
                 node_id: NodeID,
                 base: Base,
                 aligned_to: Optional[NodeID] = None,
                 column_id: Optional[ColumnID] = None,
                 block_id: Optional[BlockID] = None):
        self.node_id: NodeID = node_id
        self._base: Base = base
        self.aligned_to: NodeID = aligned_to
        self.column_id: ColumnID = column_id
        self.block_id: BlockID = block_id

    def get_base(self) -> str:
        """
        Returns base the node represents.

        Returns: Node base as str.
        """

        return self._base.value.decode("ASCII")

    def __eq__(self, other: 'Node'):
        return (self.node_id == other.node_id
                and self._base == other._base
                and self.aligned_to == other.aligned_to
                )

    def __str__(self):
        return \
            f"id: {self.node_id}, " \
            f"base: {self.get_base()}, " \
            f"aligned_to: {self.aligned_to}, " \
            f"column_id: {self.column_id}, " \
            f"block_id: {self.block_id}"

    def __repr__(self):
        return self.__str__()


SeqPath = NewType('SeqPath', List[NodeID])
SequenceMetadata = NewType('Metadata', Dict[str, Any])


class Sequence(object):
    """Sequence present in pangenome defined as SequenceID, SeqPath (list of list of nodes ids) and metadata.

    Args:
        seqid: Sequence ID.
        paths: List of SeqPaths as the sequence can be splitted into several paths.
        seqmetadata: Metadata connected with this sequence.

    Attributes:
        seqid (msa.SequenceID): Sequence ID.
        paths (List[SeqPath]): List of SeqPaths as the sequence can be splitted into several paths.
        seqmetadata (SequenceMetadata): Dictionary of metadata connected with this sequence.
    """

    def __init__(self, seqid: msa.SequenceID, paths: List[SeqPath], seqmetadata: SequenceMetadata):
        self.seqid: msa.SequenceID = seqid
        self.paths: List[SeqPath] = paths
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


class Compatibility(object):
    """Asymetric similiarity measure of two poagraph paths.

    Args:
        compatibility: Raw compatibility value - count of common nodes devided by length of one of the paths.
        p: Parameter to control compatibility value interpretation. Compatibility is raised to the power of P.

    Attributes:
        value (float): Compatibility value raised to the power of p.
        p (float): P parameter value.

    """

    def __init__(self, compatibility: float, p: at_params.P = at_params.P(1)):
        self.value: float = compatibility**p.value
        self.p: float = p.value

    def _check_p_equality(self, other: Union["Compatibility", Any]) -> None:
        if isinstance(other, Compatibility):
            assert self.p == other.p, 'Cannot compare compatibilities with different p values.'
        else:
            return

    def __eq__(self, other: Union["Compatibility", at_params.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value == other.value

    def __lt__(self, other: Union["Compatibility", at_params.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value < other.value

    def __le__(self, other: Union["Compatibility", at_params.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value <= other.value

    def __gt__(self, other: Union["Compatibility", at_params.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value > other.value

    def __ge__(self, other: Union["Compatibility", at_params.Stop]) -> bool:
        self._check_p_equality(other)
        return self.value >= other.value

    def __sub__(self, other: Union["Compatibility", at_params.Stop]) -> "Compatibility":
        self._check_p_equality(other)
        return Compatibility(self.value - other.value, at_params.P(self.p))

    def __str__(self) -> str:
        return f"""{self.value}"""

    def __repr__(self) -> str:
        return f"""value: {self.value}, p: {self.p}"""

    def base_value(self) -> "Compatibility":
        """Returns compatibility value without P transformation.

        Returns:
        Compatibility object with the original compatibility value and P=1."""

        return Compatibility(self.value ** (1 / self.p))


class Poagraph(object):
    """Pangenome is defined in Pangtree as poagraph which consists of connected nodes and sequences paths.

    Args:
        nodes: List of nodes.
        sequences: Dictonary of sequences IDs and corresponding sequences objects.
        datatype: Type of the data in poagraph.

    Attributes:
        nodes: List of nodes.
        sequences: Dictonary of sequences IDs and corresponding sequences objects.
        datatype: Type of the data in poagraph.

    """
    def __init__(self,
                 nodes: List[Node],
                 sequences: Dict[msa.SequenceID, Sequence],
                 datatype: Optional[DataType] = DataType.Nucleotides):
        self.nodes: List[Node] = nodes
        self.sequences: Dict[msa.SequenceID, Sequence] = sequences
        self.datatype: DataType = datatype

    def __eq__(self, other: 'Poagraph') -> bool:
        return self.nodes == other.nodes and \
            self.sequences == other.sequences and \
            self.datatype == other.datatype

    def get_compatibilities(self,
                            sequences_ids: List[msa.SequenceID],
                            consensus_path: SeqPath,
                            p: Optional[at_params.P] = at_params.P(1)) -> Dict[msa.SequenceID, Compatibility]:
        """Calculate compatibilities of sequences listed in sequences_ids to given consensus_path. Use P.

        Args:
            sequences_ids: Compatibilities of these seqeunces will be calculated.
            consensus_path: Sequences will be compared to this consensus path.
            p: Parameter P, see affinity tree algorithm parameters.

        Returns:
            Dictionary of sequences_ids and corresponding compatibility to given consensus_path.

        Raises:
            KeyError: If there is no sequence with given ID.
        """

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

    def get_sequences_weights(self, sequences_ids: List[msa.SequenceID]) -> \
            Dict[msa.SequenceID, int]:
        """Calculates and normalizes sequences weights inside given group of sequences.

        Args:
            sequences_ids: IDs of the sequences that the weights must be calculated and normalized.

        Returns:
            Dictionary of sequences IDs and corresponding weights.
        """

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

    def get_sequence_nodes_count(self, seq_id: msa.SequenceID) -> int:
        """Returns length of indicated sequence (sum of lengths of all sequence paths).

        Args:
            seq_id: ID of the sequence.

        Returns:
            Sum of lengths of all sequence paths.
        """

        if seq_id not in self.sequences:
            raise Exception("No sequence with given ID in poagraph.")
        return sum([len(path) for path in self.sequences[seq_id].paths])

    def get_sequences_ids(self) -> List[msa.SequenceID]:
        """Returns _sequences of all genome _sequences in Poagraph.

        Returns:
            List of sequences IDs.
        """

        return [*self.sequences.keys()]

    @staticmethod
    def complement_metadata_for_sequences_absent_in_metadata_provided(poagraph: 'Poagraph',
                                                                      metadata: msa.MetadataCSV) -> None:
        """Complements metadata in given poagraph using provided metadata.

        Args:
            poagraph: The poagraph to be changed.
            metadata: Metadata to be inserted.
        """

        headers = metadata.get_metadata_keys()
        for seq in poagraph.sequences.values():
            for h in headers:
                if h not in seq.seqmetadata:
                    seq.seqmetadata[h] = None
