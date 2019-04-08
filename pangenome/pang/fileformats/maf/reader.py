from pangenome.pang.fileformats.maf.DAGMaf import DAGMaf
from pangenome.pang.tools import loggingtools

global_logger = loggingtools.get_global_logger()


def maf_to_dagmaf(maf_content) -> DAGMaf:
    global_logger.info("Converting MAF to DAG...")
    dagmaf = DAGMaf(maf_content)
    return dagmaf
