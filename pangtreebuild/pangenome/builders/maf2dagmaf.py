from pangtreebuild.pangenome.DAGMaf import DAGMaf, DAGMafNode
from pangtreebuild.pangenome.parameters import multialignment
from pangtreebuild.mafgraph.sorter import sort_mafblocks
from pangtreebuild.tools import logprocess

global_logger = logprocess.get_global_logger()


def get_dagmaf(maf: multialignment.Maf) -> DAGMaf:
    sorted_blocks = sort_mafblocks(maf.filecontent)
    dagmafnodes = [
        DAGMafNode(block_id=b.id,
                   alignment=b.alignment,
                   orient=b.orient,
                   order=b.order(),
                   out_edges=b.out_edges)
        for b in sorted_blocks
    ]
    return DAGMaf(dagmafnodes)
