from ..context import maf_to_dagmaf
from ..context import Pangraph
from ..context import PangraphBuilderFromDAG


class TestsHelper:
    @staticmethod
    def setup_pangraph_from_maf(maf_path, metadata, fasta_source):
        dagmaf = maf_to_dagmaf(maf_path)
        pangraph = Pangraph()
        builder = PangraphBuilderFromDAG(metadata, fasta_source)
        builder.build(dagmaf, pangraph)
        return pangraph
