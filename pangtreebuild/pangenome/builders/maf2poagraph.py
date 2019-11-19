from typing import Optional, List, Tuple, Dict

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from pangtreebuild.tools import logprocess
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.pangenome import graph

_ParsedMaf = List[Optional[MultipleSeqAlignment]]

global_logger = logprocess.get_global_logger()
detailed_logger = logprocess.get_logger("details")


def get_poagraph(maf: msa.Maf, metadata: Optional[msa.MetadataCSV]) -> \
        Tuple[List[graph.Node], Dict[msa.SequenceID, graph.Sequence]]:
    """Get poagraph elements from MAF.

    Args:
        maf: Multialignment file in MAF format.
        metadata: MetadataCSV.

    Returns:
        Tuple of poagraph elements.
    """

    alignment = [*AlignIO.parse(maf.filecontent, "maf")]
    nodes, sequences = _init_poagraph(alignment, metadata)

    current_node_id = graph.NodeID(-1)
    column_id = graph.ColumnID(-1)
    for block_id, block in enumerate(alignment):
        global_logger.info(f"Processing block {block_id}...")
        block_width = len(block[0].seq)

        for col in range(block_width):
            column_id += 1
            sequence_id_to_nucleotide = {msa.SequenceID(seq.id): seq[col]
                                         for seq in block}
            nodes_codes = sorted([*(
                set([nucleotide
                     for nucleotide in sequence_id_to_nucleotide.values()])).difference({'-'})])
            column_nodes_ids = [graph.NodeID(current_node_id + i + 1)
                                for i, _ in enumerate(nodes_codes)]

            for i, nucl in enumerate(nodes_codes):
                current_node_id += 1
                nodes.append(graph.Node(node_id=current_node_id,
                                        base=graph.Base(nucl),
                                        aligned_to=_get_next_aligned_node_id(graph.NodeID(i), column_nodes_ids),
                                        column_id=graph.ColumnID(column_id),
                                        block_id=graph.BlockID(block_id)
                                        )
                             )

                for seq_id, nucleotide in sequence_id_to_nucleotide.items():
                    if nucleotide == nucl:
                        sequences[seq_id] = _add_node_do_sequence(sequence=sequences[seq_id],
                                                                  node_id=current_node_id)

    return nodes, sequences


def _init_poagraph(alignment: _ParsedMaf,
                   metadata: Optional[msa.MetadataCSV]) -> \
        Tuple[List[graph.Node], Dict[msa.SequenceID, graph.Sequence]]:
    maf_sequences_ids = _get_sequences_ids(alignment)
    metadata_sequences_ids = metadata.get_all_sequences_ids() if metadata else []
    initial_sequences = {seq_id: graph.Sequence(seqid=seq_id,
                                                paths=[],
                                                seqmetadata=metadata.get_sequence_metadata(seq_id)
                                                if metadata else {})
                         for seq_id in set(maf_sequences_ids + metadata_sequences_ids)}

    return [], initial_sequences


def _get_sequences_ids(alignment: _ParsedMaf) -> List[msa.SequenceID]:
    return list({msa.SequenceID(seq.id)
                 for block in alignment for seq in block})


def _get_next_aligned_node_id(current_column_i: graph.NodeID,
                              column_nodes_ids: List[graph.NodeID]) -> \
        Optional[graph.NodeID]:
    if len(column_nodes_ids) > 1:
        return column_nodes_ids[(current_column_i + 1) % len(column_nodes_ids)]
    return None


def _add_node_do_sequence(sequence: graph.Sequence, node_id: graph.NodeID) -> \
        graph.Sequence:
    if sequence.paths:
        a = graph.SeqPath([node_id])
        updated_path = graph.SeqPath(sequence.paths[-1] + a)
        newpaths = [sequence.paths[:-1] + updated_path]
    else:
        newpaths = [graph.SeqPath([node_id])]
    return graph.Sequence(sequence.seqid, newpaths, sequence.seqmetadata)
