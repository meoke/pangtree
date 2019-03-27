from typing import List

from pangraph.DataType import DataType
from pangraph.custom_types import NodeID


class NodePO:
    def __init__(self,
                 base: bytes,
                 aligned_to: NodeID,
                 in_nodes: List[NodeID],
                 sequences_ids: List[int]):
        self.base = base
        self.aligned_to = aligned_to
        self.in_nodes = in_nodes
        self.sequences_ids = sequences_ids


class SequencePO:
    def __init__(self, name: str, nodes_count: int, weight: int, consensus_id: int, start_node_id: NodeID):
        self.name = name
        self.nodes_count = nodes_count
        self.weight = weight
        self.consensus_id = consensus_id
        self.start_node_id = start_node_id


class PangraphToPO:
    def __init__(self):
        self.po_nodes: List[NodePO] = None
        self.po_sequences: List[SequencePO] = None
        self.po_lines: List[str] = None

    def get_po_file_content(self, po_nodes: List[NodePO], po_sequences: List[SequencePO], datatype: DataType) -> str:
        self.po_nodes = po_nodes
        self.po_sequences = po_sequences
        self.po_lines = [""] * self.get_po_file_lines_count()

        last_position = self._write_introduction()
        last_position = self._write_sequences_info(start_at=last_position+1)
        _ = self._write_nodes_info(start_at=last_position+1, datatype=datatype)

        return "\n".join(self.po_lines)

    def get_po_file_lines_count(self) -> int:
        const_lines_count = 5
        source_lines_count = 2 * len(self.po_sequences)
        nodes_lines_count = len(self.po_nodes)
        return const_lines_count + source_lines_count + nodes_lines_count

    def _write_introduction(self) -> int:
        self.po_lines[0] = "VERSION=pangenome"
        self.po_lines[1] = "NAME=pangenome"
        self.po_lines[2] = "TITLE=pangenome"
        self.po_lines[3] = f"LENGTH={len(self.po_nodes)}"
        self.po_lines[4] = f"SOURCECOUNT={len(self.po_sequences)}"
        return 4

    def _write_sequences_info(self, start_at: int) -> int:
        if not self.po_sequences:
            raise Exception("No sequences info to write in PO file.")
        i = -2
        for sequence in self.po_sequences:
            i += 2
            self.po_lines[start_at + i] = f"SOURCENAME={sequence.name}"
            self.po_lines[start_at + i + 1] = ("SOURCEINFO=" +
                                               " ".join([f"{sequence.nodes_count}",
                                                         f"{sequence.start_node_id}",
                                                         f"{sequence.weight}",
                                                         f"{sequence.consensus_id}",
                                                         f"{sequence.name}"]))
        return start_at + i + 1

    def _write_nodes_info(self, start_at: int, datatype: DataType) -> int:
        if not self.po_nodes:
            raise Exception("No nodes info to write in PO file.")
        i = 0

        if datatype.value == DataType.Proteins.value:
            _get_node_code = PangraphToPO._get_protein_node_code
        elif datatype.value == DataType.Nucleotides.value:
            _get_node_code = PangraphToPO._get_nucleotides_node_code
        else:
            raise Exception("Unknown data type. Cannot create PO file.")

        for node in self.po_nodes:
            self.po_lines[start_at + i] = "".join([_get_node_code(node.base),
                                                   ":",
                                                   self._get_in_nodes_info(node.in_nodes),
                                                   self._get_sources_info(node.sequences_ids),
                                                   self._get_aligned_to_info(node.aligned_to)])
            i += 1
        return start_at + i

    @staticmethod
    def _get_protein_node_code(base: bytes) -> str:
        return base.decode("ASCII").upper()

    @staticmethod
    def _get_nucleotides_node_code(base: bytes) -> str:
        return base.decode("ASCII").lower()

    def _get_in_nodes_info(self, in_nodes: List[NodeID]) -> str:
        return "".join([f'L{i}' for i in in_nodes])

    def _get_sources_info(self, sources_ids: List[int]) -> str:
        return "".join([f'S{i}' for i in sources_ids])

    def _get_aligned_to_info(self, aligned_to: NodeID) -> str:
        return f"A{aligned_to}" if aligned_to is not None else ""
