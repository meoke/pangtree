from pangtreebuild.mafgraph.sorter import sort_mafblocks
from pangtreebuild.pangenome import DAGMaf
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.tools import logprocess

global_logger = logprocess.get_global_logger()


def get_dagmaf(maf: msa.Maf) -> DAGMaf.DAGMaf:
    """Converts MAF to DagMaf.

    Args:
        maf: MAF to be converted.

    Returns:
         DagMaf built from the MAF.
    """
    sorted_blocks = sort_mafblocks(maf.filecontent)
    dagmafnodes = [
        DAGMaf.DAGMafNode(block_id=b.id,
                          alignment=b.alignment,
                          orient=b.orient,
                          order=b.order(),
                          out_edges=b.out_edges)
        for b in sorted_blocks
    ]
    return DAGMaf.DAGMaf(dagmafnodes)
