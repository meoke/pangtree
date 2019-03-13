from fileformats.maf.DAGMaf import DAGMaf


def maf_to_dagmaf(maf_content) -> DAGMaf:
    dagmaf = DAGMaf(maf_content)
    return dagmaf
