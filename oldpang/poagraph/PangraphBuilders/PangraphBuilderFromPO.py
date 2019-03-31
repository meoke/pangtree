from collections import namedtuple
from io import StringIO
from typing import Dict, Iterator, List, Tuple

from data.Node import Node
from data.PangraphBuilders.PangraphBuilderBase import PangraphBuilderBase
from metadata.MultialignmentMetadata import MultialignmentMetadata
from data.custom_types import NodeID, SequenceID, Base, make_base
from tools import loggingtools

global_logger = loggingtools.get_global_logger()
detailed_logger = loggingtools.get_logger("details")

POSequenceInfo = namedtuple('POSequenceInfo', ['name',
                                               'nodes_count',
                                               'start_node',
                                               'weight',
                                               'consensus_id',
                                               'additional_info'])


class PangraphBuilderFromPO(PangraphBuilderBase):
    def __init__(self, genomes_info: MultialignmentMetadata, missing_base_symbol):
        super().__init__(genomes_info)
        self.pangraph = None
        self.missing_nucleotide_symbol: Base = make_base(missing_base_symbol)
        self.sequences_info: Dict[int, POSequenceInfo] = {}

    def build(self, po_content: StringIO, pangraph):
        self.pangraph = pangraph
        for i in range(3):
            next(po_content)

        nodes_count = int(PangraphBuilderFromPO.extract_line_value(next(po_content)))
        paths_count = int(PangraphBuilderFromPO.extract_line_value(next(po_content)))

        self._init_pangraph(nodes_count)
        self._build_sequences_info(po_content, paths_count)
        self._build_pangraph_paths_and_nodes(po_content, nodes_count)

    def _init_pangraph(self, nodes_count: int) -> None:
        self.pangraph.paths = {SequenceID(seq_id): [[]] for seq_id in self.sequences_ids}
        self.pangraph.nodes = [None] * nodes_count

    @staticmethod
    def extract_line_value(line: str) -> str:
        splitted_line = line.split('=')
        if len(splitted_line) < 2:
            raise Exception("Incorrect line format. Nothing after \"=\" found.")
        return splitted_line[1].strip()

    @staticmethod
    def get_sequences_names(multialignment_file_content: str) -> List[str]:
        pocontent = StringIO(multialignment_file_content)
        po_lines = pocontent.readlines()
        sourcecount = int(PangraphBuilderFromPO.extract_line_value(po_lines[4]))
        sequences_names = []
        for i in range(5, 5+sourcecount*2, 2):
            sequence_name = PangraphBuilderFromPO.extract_line_value(po_lines[i])
            sequences_names.append(sequence_name)
        return sequences_names

    def _build_sequences_info(self, po_content: Iterator, paths_count: int) -> None:
        for i in range(paths_count):
            path_name = PangraphBuilderFromPO.extract_line_value(next(po_content))
            detailed_path_info = next(po_content)
            detailed_info = PangraphBuilderFromPO.extract_line_value(detailed_path_info).split(' ')
            if len(detailed_info) < 5:
                raise Exception(f"Expeceted SOURCEINFO=[5 parameters]. Got {detailed_path_info} instead.")
            self.sequences_info[i] = POSequenceInfo(name=path_name,
                                                    nodes_count=detailed_info[0],
                                                    start_node=detailed_info[1],
                                                    weight=detailed_info[2],
                                                    consensus_id=detailed_info[3],
                                                    additional_info="".join(detailed_info[4:]))

    def _build_pangraph_paths_and_nodes(self, po_content: Iterator, nodes_count: int) -> None:
        for i in range(nodes_count):
            node_line = next(po_content)
            base = make_base(node_line[0].upper())
            in_nodes, po_sequences_ids, aligned_to = PangraphBuilderFromPO._extract_node_parameters(node_line)
            sequences_ids = [self.sequences_info[po_sequences_id].name for po_sequences_id in po_sequences_ids]
            self.add_node_to_pangraph(NodeID(i), base, NodeID(aligned_to), sequences_ids)

    def add_node_to_pangraph(self, node_id: NodeID, base: Base, aligned_to: NodeID, seq_ids: List[SequenceID]) -> None:
        self.pangraph.nodes[node_id] = Node(node_id=node_id,
                                            base=base,
                                            aligned_to=aligned_to
                                            )
        for seq_id in seq_ids:
            self.pangraph.paths[seq_id][0].append(node_id)

    @staticmethod
    def _extract_node_parameters(line: str) -> Tuple[List[int], List[int], int]:
        line = line[2:]
        line_iter = iter(line)

        node_parameters = {'L': [], 'S': [], 'A': [], '': []}
        current_node_parameter = next(line_iter)
        number_start = 1
        number_end = 1
        for i, c in enumerate(line_iter):
            if c == 'L' or c == 'S' or c == 'A':
                node_parameters[current_node_parameter].append(int(line[number_start: number_end]))
                current_node_parameter = c
                number_start = i + 2
                number_end = i + 2
            else:
                number_end += 1
        node_parameters[current_node_parameter].append(int(line[number_start: number_end]))
        aligned_to = node_parameters['A'][0] if node_parameters['A'] else None
        return node_parameters['L'], node_parameters['S'], aligned_to
