from pangtreebuild.datamodel.DAGMaf import DAGMaf, DAGMafNode
from pangtreebuild.datamodel.input_types import Maf
from mafgraph.sorter import sort_mafblocks
from pangtreebuild.tools import logprocess

global_logger = logprocess.get_global_logger()


def get_dagmaf(maf: Maf) -> DAGMaf:
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
