import os
from bisect import bisect_left
from collections import namedtuple
from pathlib import Path
from typing import List, Dict, Union

from datamodel.Node import NodeID
from datamodel.Poagraph import Poagraph
from datamodel.Sequence import SequenceID, SequencePath
from output.PangenomePO import NodePO, SequencePO
from tools import pathtools
import output.PangenomePO as PangenomePO
import subprocess


class NoConsensusError(Exception):
    pass


def get_consensus(poagraph: Poagraph,
                  sequences_ids: List[SequenceID],
                  output_dir: Path,
                  job_name: str,
                  blosum_path: Path,
                  specific_consensus_id: int) -> Dict[int, SequencePath]:
    poa_input_path = pathtools.get_child_path(output_dir, f"{job_name}_in_pangenome.po")
    poa_output_path = pathtools.get_child_path(output_dir, f"{job_name}_out_pangenome.po")

    s = PangraphPOTranslator(poagraph, sequences_ids)
    poa_input_content = s.get_input_po_content()
    with open(poa_input_path, 'w') as poa_input:
        poa_input.write(poa_input_content)
    b_resolved = blosum_path.resolve()
    call(po_file_path=poa_input_path,
         hb_file_path=poa_output_path,
         blosum_path=blosum_path.resolve(),
         hbmin=0.6)
    with open(poa_output_path) as poa_output:
        poa_output_lines = poa_output.readlines()
    consensus_paths = s.read_consensus_paths(poa_output_lines, [specific_consensus_id])
    return consensus_paths
    # return {specific_consensus_id: SequencePath([NodeID(0)])}


def call(po_file_path: Path, hb_file_path: Path, blosum_path: Path, hbmin: float) -> None:
    poa_path = pathtools.get_child_path(Path(os.path.abspath(__file__)), '../../../bin/poa').resolve()
    # detailed_logger.info(f"Run poa! Input: {po_file_path} Output: {hb_file_path}...")
    command = f"{poa_path} -read_msa {po_file_path} -hb -po {hb_file_path} {blosum_path} -hbmin {hbmin}"
    # poa_result = subprocess.run(command, stderr=subprocess.PIPE, shell=True)
    poa_result = subprocess.run(command, stderr=subprocess.PIPE, shell=True)
    poa_str_output = poa_result.stderr.decode("ASCII")
    # detailed_logger.info(f"Poa output: {poa_str_output}")



class PangraphPOTranslator:
    def __init__(self, poagraph: Poagraph, sequences_ids: List[SequenceID]):
        self.poagraph: Poagraph = poagraph
        self.sequences_ids: List[SequenceID] = sequences_ids
        self.new_to_old: Dict[NodeID, NodeID] = None
        self.old_to_new: Dict[NodeID, NodeID] = None
        self.seq_old_to_new: Dict[SequenceID, int] = None
        self.seq_new_to_old: Dict[int, SequenceID] = None

    def get_input_po_content(self) -> str:
        paths_to_keep = [path
                         for seq_id in self.sequences_ids
                         for path in self.poagraph.sequences[seq_id].paths]
        nodes_ids_to_keep = list(set([node_id
                                      for path in paths_to_keep
                                      for node_id in path]))
        sorted_nodes_ids_to_keep = sorted(nodes_ids_to_keep)
        self.old_to_new = {node_id: i for i, node_id in enumerate(sorted_nodes_ids_to_keep)}
        self.new_to_old = {new_node_id: old_node_id for old_node_id, new_node_id in self.old_to_new.items()}
        self.seq_old_to_new = {seq_id: i for i, seq_id in enumerate(self.sequences_ids)}
        self.seq_new_to_old = {i: seq_id for seq_id, i in self.seq_old_to_new.items()}

        po_nodes = [NodePO(base=self.poagraph.nodes[self.new_to_old[new_node_id]].base.value,
                           aligned_to=self._get_aligned_node(self.new_to_old[new_node_id], sorted_nodes_ids_to_keep),
                           in_nodes=set(),
                           sequences_ids=[]
                           )
                    for new_node_id in range(len(nodes_ids_to_keep))]

        sequences_weight: Dict[SequenceID, int] = self.poagraph.get_sequences_weights(self.sequences_ids)
        po_sequences = [SequencePO(name=self.seq_new_to_old[new_seq_id],
                                   nodes_count=self.poagraph.get_sequence_nodes_count(self.seq_new_to_old[new_seq_id]),
                                   weight=sequences_weight[self.seq_new_to_old[new_seq_id]],
                                   consensus_id=-1,
                                   start_node_id=self.old_to_new[
                                       self.poagraph.sequences[self.seq_new_to_old[new_seq_id]].paths[0][0]]
                                   )
                        for new_seq_id in range(len(self.sequences_ids))]

        for seq_id in self.sequences_ids:
            new_seq_id = self.seq_old_to_new[seq_id]
            for path in self.poagraph.sequences[seq_id].paths:
                for i, node_id in enumerate(path):
                    new_node_id = self.old_to_new[node_id]
                    po_nodes[new_node_id].sequences_ids.append(new_seq_id)
                    if i > 0:
                        new_in_node_id = self.old_to_new[path[i-1]]
                        po_nodes[new_node_id].in_nodes.add(new_in_node_id)

        for node in po_nodes:
            node.in_nodes = sorted(node.in_nodes)

        return PangenomePO.poagraph_elements_to_PangenomePO(po_nodes, po_sequences, self.poagraph.datatype)

    def _get_aligned_node(self, old_node_id: NodeID, sorted_nodes_ids_to_keep: List[NodeID]) -> Union[NodeID, None]:
        aligned_to = self.poagraph.nodes[old_node_id].aligned_to
        if aligned_to is None:
            return None
        while aligned_to != old_node_id:
            if self._is_in(sorted_list=sorted_nodes_ids_to_keep, x=aligned_to):
                return self.old_to_new[aligned_to]
            aligned_to = self.poagraph.nodes[aligned_to].aligned_to
        return None

    def _is_in(self, sorted_list, x):
        i = bisect_left(sorted_list, x)
        if i != len(sorted_list) and sorted_list[i] == x:
            return True
        return False

    @staticmethod
    def _extract_line_value(line: str) -> str:
        return line.split('=')[1].strip()

    def read_consensus_paths(self, po_lines: List[str], specific_consensus_id) -> Dict[int, SequencePath]:
        po_lines_iterator = iter(po_lines)

        for i in range(3):
            next(po_lines_iterator)

        nodes_count = int(PangraphPOTranslator._extract_line_value(next(po_lines_iterator)))
        paths_count = int(PangraphPOTranslator._extract_line_value(next(po_lines_iterator)))

        ConsInfo = namedtuple('ConsInfo', ['fullname', 'po_consensus_id', 'length'])

        consensuses_in_po_lines = dict()
        s = -1
        for i in range(paths_count):
            s += 1
            path_name = PangraphPOTranslator._extract_line_value(next(po_lines_iterator))
            if len(path_name) > 7 and path_name[:7] == "CONSENS":
                detailed_consens_info = next(po_lines_iterator)
                detailed_info = self._extract_line_value(detailed_consens_info).split(' ')
                consens_nodes_count = int(detailed_info[0])
                consensuses_in_po_lines[int(path_name[7:])] = ConsInfo(fullname=path_name,
                                                                       po_consensus_id=f"S{str(s)}",
                                                                       length=consens_nodes_count)
            else:
                _ = next(po_lines_iterator)
                continue

        if not consensuses_in_po_lines:
            raise NoConsensusError("No consensus found in this poagraph.")

        consensuses_to_return = {c: [None] * c_info.length
                                 for c, c_info
                                 in consensuses_in_po_lines.items()
                                 if c in specific_consensus_id}
        consensuses_paths_node_counter = {c: 0 for c in consensuses_to_return.keys()}

        first_node_line_idx = 5 + paths_count * 2
        old_node_id = 0
        new_node_id = -1

        for line_idx in range(first_node_line_idx, len(po_lines)):
            new_node_id += 1
            for c_id in consensuses_to_return.keys():
                if consensuses_in_po_lines[c_id].po_consensus_id in po_lines[line_idx]:
                    consensuses_to_return[c_id][consensuses_paths_node_counter[c_id]] = self.new_to_old[new_node_id]
                    consensuses_paths_node_counter[c_id] += 1
        return consensuses_to_return