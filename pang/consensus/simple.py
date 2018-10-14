from graph.Pangraph import Pangraph
import po.writer as powriter
import po.reader as poreader
import userio.pathtools as pathtools
from pathlib import Path
import consensus.poa as poa
from metadata.MultialignmentMetadata import MultialignmentMetadata


def run(outputdir: Path, pangraph: Pangraph, hbmin: float, genomes_info: MultialignmentMetadata) -> Pangraph:
    poa_input_path = pathtools.get_child_file_path(outputdir, "in_pangenome.po")
    poa_output_path = pathtools.get_child_file_path(outputdir, "out_pangenome.po")
    powriter.save(pangraph, poa_input_path, genomes_info=genomes_info)

    poa.run(input_path=poa_input_path, output_path=poa_output_path, hbmin=hbmin)

    pangraph = poreader.read(poa_output_path)
    return pangraph

#
# def get_top_consensus(poagraph, sources_IDs, hbmin):
#     po_file_path, nodes_map = po_writer.save_as_po(poagraph, sources_IDs)
#     hb_file_path = t.change_file_extension(po_file_path, '.hb')
#     run(['../bin/poa', '-read_msa', po_file_path, '-hb', '-po', hb_file_path, '../bin/blosum80.mat', '-v', '-hbmin',
#          str(hbmin)])
#     try:
#         consensus0, consensus_nodes = po_reader.read_consensus(hb_file_path, consensusID=0)
#     except NoConsensusFound:
#         raise NoConsensusFound()
#
#     consensus_actual_nodes = np.zeros(shape=len(poagraph.nodes), dtype=np.bool)
#     for i, val in enumerate(consensus_nodes):
#         orig_ID = nodes_map['orig_ID'][nodes_map['temp_ID'] == i][0]
#         consensus_actual_nodes[orig_ID] = val
#
#     return consensus0, consensus_actual_nodes