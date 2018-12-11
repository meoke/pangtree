from graph.Pangraph import Pangraph
import po.writer as powriter
import po.reader as poreader
import userio.pathtools as pathtools
from pathlib import Path
import consensus_algorithm.poa as poa
from metadata.MultialignmentMetadata import MultialignmentMetadata


def run(outputdir: Path, pangraph: Pangraph, hbmin: float, genomes_info: MultialignmentMetadata, filename_prefix: str) -> Pangraph:
    poa_input_path = pathtools.get_child_file_path(outputdir, f"{filename_prefix}_in_pangenome.po")
    poa_output_path = pathtools.get_child_file_path(outputdir, f"{filename_prefix}_out_pangenome.po")
    powriter.save(pangraph, poa_input_path, genomes_info=genomes_info)

    poa.run(input_path=poa_input_path, output_path=poa_output_path, hbmin=hbmin)

    pangraph = poreader.read(poa_output_path)
    return pangraph
