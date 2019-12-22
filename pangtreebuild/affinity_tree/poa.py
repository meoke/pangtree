"""Functions to call poa software."""

from bisect import bisect_left
import os
from pathlib import Path
import subprocess
from typing import List, Dict, Union, Optional

from pangtreebuild.affinity_tree import parameters
from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.serialization import po
from pangtreebuild.tools import pathtools
from pangtreebuild.tools import logprocess


detailed_logger = logprocess.get_logger('details')
global_logger = logprocess.get_global_logger()


class NoConsensusError(Exception):
    """If poa software cannot find any consensus path in the Poagraph"""

    pass


class ConsInfo(object):
    """Consensus information based on poa result.

    Args:
        fullname: Consensus path name. Eg. "CONSENS0"
        po_consensus_id: Consensus ID used in PO file to indicate it is
            present in a node. Eg. "S0"
        assigned_sequences_ids: IDs of _sequences assigned to this consensus.
        path: List of nodes of this consensus.

    Attributes:
        fullname (str): Consensus name.
        po_consensus_id (str): Consensus ID used in PO file to indicate
            it is present in a node. Eg. "S0"
        assigned_sequences_ids (List[msa.SequenceID]): IDs of sequences
            assigned to this consensus.
        path (graph.SeqPath): List of nodes of this consensus.
    """

    def __init__(self,
                 fullname: str,
                 po_consensus_id: Optional[str] = None,
                 assigned_sequences_ids: Optional[List[msa.SequenceID]] = None,
                 path: Optional[graph.SeqPath] = None):
        self.fullname: str = fullname
        self.po_consensus_id: str = po_consensus_id
        self.assigned_sequences_ids: List[msa.SequenceID] = assigned_sequences_ids
        self.path: graph.SeqPath = path


def get_consensuses(poagraph: graph.Poagraph,
                    sequences_ids: List[msa.SequenceID],
                    output_dir: Path,
                    job_name: str,
                    blosum_path: Path,
                    hbmin: parameters.Hbmin,
                    specific_consensuses_id: Optional[List[int]] = None) -> \
                        Dict[int, ConsInfo]:
    """Calls poa software on given Poagraph to get consensus paths.

    Args:
        poagraph: Poagraph used as input to poa software. It may be cropped by
            using sequences_ids argument.
        sequences_ids: IDs of the _sequences that should be kept in poagraph
            being input to poa.
        output_dir: Full path to the directory used by poa software as
            temporary storage place.
        job_name: Name of the task used to label produced file names.
        blosum_path: Full path to the Blosum file used as poa's input.
        hbmin: Hbmin value used used as poa's input.
        specific_consensuses_id: Poa returns consensuses numbered by: 0, 1...
            By this argument it can be specified which should be returned.

    Returns:
        Dictionary of consensus numbers and corresponding information as
            ConsInfo object.

    Raises:
        NoConsensusError: If no consensus was found for given Poagraph and
            set of selected _sequences.
    """
    poa_input_path = pathtools.get_child_path(output_dir,
                                              f"{job_name}_in_pangenome.po")
    poa_output_path = pathtools.get_child_path(output_dir,
                                               f"{job_name}_out_pangenome.po")

    s = _PoagraphPOTranslator(poagraph, sequences_ids)
    poa_input_content = s.get_input_po_content()
    with open(poa_input_path, 'w') as poa_input:
        poa_input.write(poa_input_content)
    _call_poa(po_file_path=poa_input_path,
              hb_file_path=poa_output_path,
              blosum_path=blosum_path.resolve(),
              hbmin=hbmin.value)
    with open(poa_output_path) as poa_output:
        poa_output_lines = poa_output.readlines()
    os.remove(poa_input_path)
    os.remove(poa_output_path)
    consensus_paths = s.read_consensus_paths(poa_output_lines,
                                             specific_consensuses_id)
    return consensus_paths


def _call_poa(po_file_path: Path,
              hb_file_path: Path,
              blosum_path: Path,
              hbmin: float) -> None:
    """Calls poa software.

    Args:
        po_file_path: Path to the PO file containing Poagraph.
        hb_file_path: Path to the PO file where the poa result will be saved.
        blosum_path: Path to the Blosum file.
        hbmin: Hbmin value used as poa input.

    Returns:
        Nothing, as the result is stored in hb_file_path.
    """
    affinity_tree_dir = Path(__file__).parent
    poa_path = pathtools.get_child_path(affinity_tree_dir, "bin/poa")
    detailed_logger.info(f"Run poa! Input: {po_file_path} Output: {hb_file_path}...")
    command = f"{poa_path} -read_msa {po_file_path} -hb -po {hb_file_path} {blosum_path} -hbmin {hbmin}"
    poa_result = subprocess.run(command, stderr=subprocess.PIPE, shell=True)
    poa_str_output = poa_result.stderr.decode("ASCII")
    detailed_logger.info(f"Poa output: {poa_str_output}")


class _PoagraphPOTranslator:
    """Converts Poagraph to PO file and back poa result PO file to Poagraph.

    Args:
        poagraph: Full poagraph where the consensuses are searched.
        sequences_ids: Sequences to be kept in the final PO file
            which stores Poagraph.

    Attributes:
        poagraph (poagraph.Poagraph): The orginal poagraph.
        sequences_ids (List[msa.SequenceID]): List of _sequences IDs to be
            included in poa input.
        new_to_old: Dict[poagraph.NodeID, poagraph.NodeID]: Mapping of
            temporary poagraph nodes IDs to the original ones.
        old_to_new: Dict[poagraph.NodeID, poagraph.NodeID] = Mapping of the
            original poagraph nodes IDs to the temporary ones.
        seq_old_to_new: Dict[msa.SequenceID, int] = Mapping of the original
            sequences IDs to temporary integer ones.
        seq_new_to_old: Dict[int, msa.SequenceID] = Mapping of the temporary
            int IDs of _sequences to the original ones.
    """
    def __init__(self,
                 poagraph: graph.Poagraph,
                 sequences_ids: List[msa.SequenceID]):
        self.poagraph: graph.Poagraph = poagraph
        self.sequences_ids: List[msa.SequenceID] = sequences_ids
        self.new_to_old: Dict[graph.NodeID, graph.NodeID] = None
        self.old_to_new: Dict[graph.NodeID, graph.NodeID] = None
        self.seq_old_to_new: Dict[msa.SequenceID, int] = None
        self.seq_new_to_old: Dict[int, msa.SequenceID] = None

    def get_input_po_content(self) -> str:
        """Convert poagraph to PO file content."""

        paths_to_keep = [path
                         for seq_id in self.sequences_ids
                         for path in self.poagraph.sequences[seq_id].paths]
        nodes_ids_to_keep = list(set([node_id
                                      for path in paths_to_keep
                                      for node_id in path]))
        sorted_nodes_ids_to_keep = sorted(nodes_ids_to_keep)
        self.old_to_new = {node_id: i
                           for i, node_id in enumerate(sorted_nodes_ids_to_keep)}
        self.new_to_old = {new_node_id: old_node_id
                           for old_node_id, new_node_id in self.old_to_new.items()}
        self.seq_old_to_new = {seq_id: i
                               for i, seq_id in enumerate(self.sequences_ids)}
        self.seq_new_to_old = {i: seq_id
                               for seq_id, i in self.seq_old_to_new.items()}

        po_nodes = [po.NodePO(base=self.poagraph.nodes[self.new_to_old[new_node_id]]._base.value,
                              aligned_to=self._get_aligned_node(self.new_to_old[new_node_id],
                                                                sorted_nodes_ids_to_keep),
                              in_nodes=set(),
                              sequences_ids=[])
                    for new_node_id in range(len(nodes_ids_to_keep))]

        sequences_weight: Dict[msa.SequenceID, int] = self.poagraph.get_sequences_weights(self.sequences_ids)
        po_sequences = [po.SequencePO(name=self.seq_new_to_old[new_seq_id],
                                      nodes_count=self.poagraph.get_sequence_nodes_count(self.seq_new_to_old[new_seq_id]),
                                      weight=sequences_weight[self.seq_new_to_old[new_seq_id]],
                                      consensus_id=-1,
                                      start_node_id=self.old_to_new[
                                          self.poagraph.sequences[self.seq_new_to_old[new_seq_id]].paths[0][0]])
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

        return po.poagraph_elements_to_PangenomePO(po_nodes,
                                                   po_sequences,
                                                   self.poagraph.datatype)

    def _get_aligned_node(self,
                          old_node_id: graph.NodeID,
                          sorted_nodes_ids_to_keep: List[graph.NodeID]) -> \
            Union[graph.NodeID, None]:
        """Get aligned node ID or None if it there is no aligned node.

        Args:
            old_node_id: Node ID in the original poagraph.
            sorted_nodes_ids_to_keep: Node IDs that are kept in the
                cropped Poagraph.

        Returns:
            Aligned node ID if it is present or None if not.
        """
        aligned_to = self.poagraph.nodes[old_node_id].aligned_to
        if aligned_to is None:
            return None
        while aligned_to != old_node_id:
            if self._is_in(sorted_list=sorted_nodes_ids_to_keep, x=aligned_to):
                return self.old_to_new[aligned_to]
            aligned_to = self.poagraph.nodes[aligned_to].aligned_to
        return None

    def _is_in(self, sorted_list: List[int], x: int) -> bool:
        """Checks if x is in sorted list.

        Args:
            sorted_list: List of sorted numerical values.
            x: Numerical value that is searched in sorted_list.

        Returns:
            Information whether x is in sorted_list or not.
        """

        i = bisect_left(sorted_list, x)
        if i != len(sorted_list) and sorted_list[i] == x:
            return True
        return False

    @staticmethod
    def _extract_line_value(line: str) -> str:
        """Specific to PO format. Extracts value assigned to variable in line.

        Args:
            Line of the PO file in format [...]=[...]

        Returns:
            Assigned value as str."""
        return line.split('=')[1].strip()

    def read_consensus_paths(self,
                             po_lines: List[str],
                             specific_consensuses_ids: Optional[List[int]] = None) -> Dict[int, ConsInfo]:
        """For given lines of PO file reads consensus information.

        Args:
            po_lines: Lines of the PO file that is poa result.
            specific_consensuses_ids: Poa generates as many consensuses as is possible but the number of read
                                      consensuses can be modifed by specifing IDs of consensuses to return.

        Returns:
            Dictionary of consensuses IDs and corresponding Consensus Info objects.
        """
        po_lines_iterator = iter(po_lines)

        for i in range(3):
            next(po_lines_iterator)

        nodes_count = int(_PoagraphPOTranslator._extract_line_value(next(po_lines_iterator)))
        paths_count = int(_PoagraphPOTranslator._extract_line_value(next(po_lines_iterator)))

        consensuses_in_po_lines: Dict[int, ConsInfo] = dict()
        current_consensus_id = -1
        for i in range(paths_count):
            path_name = _PoagraphPOTranslator._extract_line_value(next(po_lines_iterator))
            if len(path_name) > 7 and path_name[:7] == "CONSENS":
                current_consensus_id += 1
                if specific_consensuses_ids is None or current_consensus_id in specific_consensuses_ids:
                    detailed_consens_info_line = next(po_lines_iterator)
                    detailed_consens_info = self._extract_line_value(detailed_consens_info_line).split(' ')
                    consens_nodes_count = int(detailed_consens_info[0])
                    cons_id = int(path_name[7:])
                    if cons_id in consensuses_in_po_lines:
                        consensuses_in_po_lines[cons_id].fullname = path_name
                        consensuses_in_po_lines[cons_id].po_consensus_id = f"S{str(i)}"
                        consensuses_in_po_lines[cons_id].path = [None] * consens_nodes_count
                    else:
                        consensuses_in_po_lines[cons_id] = ConsInfo(fullname=f"CONSENS{cons_id}",
                                                                    assigned_sequences_ids=[],
                                                                    po_consensus_id=f"S{str(i)}",
                                                                    path=[None] * consens_nodes_count)
            else:
                detailed_sequence_info_line = next(po_lines_iterator)
                detailed_sequence_info = self._extract_line_value(detailed_sequence_info_line).split(' ')
                assigned_consensus_id = int(detailed_sequence_info[3])
                sequence_id = self.seq_new_to_old[i]
                if assigned_consensus_id != -1:
                    if assigned_consensus_id in consensuses_in_po_lines:
                        consensuses_in_po_lines[assigned_consensus_id].assigned_sequences_ids.append(sequence_id)
                    else:
                        consensuses_in_po_lines[assigned_consensus_id] = ConsInfo(fullname=f"CONSENS{assigned_consensus_id}",
                                                                                  assigned_sequences_ids=[sequence_id])

        if not consensuses_in_po_lines:
            raise NoConsensusError("No consensus found in this poagraph.")

        consensuses_paths_node_counter = {c: 0
                                          for c in consensuses_in_po_lines.keys()}

        first_node_line_idx = 5 + paths_count * 2
        new_node_id = -1

        for line_idx in range(first_node_line_idx, len(po_lines)):
            new_node_id += 1
            for c_id in range(len(consensuses_in_po_lines)):
                if consensuses_in_po_lines[c_id].po_consensus_id in po_lines[line_idx]:
                    consensuses_in_po_lines[c_id].path[consensuses_paths_node_counter[c_id]] = self.new_to_old[new_node_id]
                    consensuses_paths_node_counter[c_id] += 1

        return consensuses_in_po_lines
