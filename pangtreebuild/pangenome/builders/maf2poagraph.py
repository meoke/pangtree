from typing import Optional, List, Tuple, Dict

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

from pangtreebuild.tools import logprocess
from pangtreebuild.pangenome.parameters import multialignment
from pangtreebuild.pangenome import poagraph

_ParsedMaf = List[Optional[MultipleSeqAlignment]]

global_logger = logprocess.get_global_logger()
detailed_logger = logprocess.get_logger("details")


def get_poagraph(maf: multialignment.Maf, metadata: Optional[multialignment.MetadataCSV]) -> Tuple[List[poagraph.Node], Dict[multialignment.SequenceID, poagraph.Sequence]]:
    alignment = [*AlignIO.parse(maf.filecontent, "maf")]
    nodes, sequences = _init_poagraph(alignment, metadata)

    current_node_id = poagraph.NodeID(-1)
    column_id = poagraph.ColumnID(-1)
    for block_id, block in enumerate(alignment):
        global_logger.info(f"Processing block {block_id}...")
        block_width = len(block[0].seq)

        for col in range(block_width):
            column_id += 1
            sequence_id_to_nucleotide = {multialignment.SequenceID(seq.id): seq[col] for seq in block}
            nodes_codes = sorted([*(
                set([nucleotide for nucleotide in sequence_id_to_nucleotide.values()])).difference({'-'})])
            column_nodes_ids = [poagraph.NodeID(current_node_id + i + 1) for i, _ in enumerate(nodes_codes)]

            for i, nucl in enumerate(nodes_codes):
                current_node_id += 1
                nodes.append(poagraph.Node(node_id=current_node_id,
                                  base=poagraph.Base(nucl),
                                  aligned_to=_get_next_aligned_node_id(poagraph.NodeID(i), column_nodes_ids),
                                  column_id=poagraph.ColumnID(column_id),
                                  block_id=poagraph.BlockID(block_id)
                                  )
                             )

                for seq_id, nucleotide in sequence_id_to_nucleotide.items():
                    if nucleotide == nucl:
                        sequences[seq_id] = _add_node_do_sequence(sequence=sequences[seq_id],
                                                                  node_id=current_node_id)

    return nodes, sequences


def _init_poagraph(alignment: _ParsedMaf,
                   metadata: Optional[multialignment.MetadataCSV]) -> Tuple[List[poagraph.Node], Dict[multialignment.SequenceID, poagraph.Sequence]]:
    maf_sequences_ids = _get_sequences_ids(alignment)
    metadata_sequences_ids = metadata.get_all_sequences_ids() if metadata else []
    initial_sequences = {seq_id: poagraph.Sequence(seqid=seq_id,
                                          paths=[],
                                          seqmetadata=metadata.get_sequence_metadata(seq_id)
                                          if metadata else {})
                         for seq_id in set(maf_sequences_ids + metadata_sequences_ids)}

    return [], initial_sequences


def _get_sequences_ids(alignment: _ParsedMaf) -> List[multialignment.SequenceID]:
    return list({multialignment.SequenceID(seq.id) for block in alignment for seq in block})


def _get_next_aligned_node_id(current_column_i: poagraph.NodeID, column_nodes_ids: List[poagraph.NodeID]) -> Optional[poagraph.NodeID]:
    if len(column_nodes_ids) > 1:
        return column_nodes_ids[(current_column_i + 1) % len(column_nodes_ids)]
    return None


def _add_node_do_sequence(sequence: poagraph.Sequence, node_id: poagraph.NodeID) -> poagraph.Sequence:
    if sequence.paths:
        a = poagraph.SeqPath([node_id])
        updated_path = poagraph.SeqPath(sequence.paths[-1] + a)
        newpaths = [sequence.paths[:-1] + updated_path]
    else:
        newpaths = [poagraph.SeqPath([node_id])]
    return poagraph.Sequence(sequence.seqid, newpaths, sequence.seqmetadata)
