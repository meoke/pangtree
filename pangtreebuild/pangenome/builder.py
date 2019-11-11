from typing import Optional, Tuple

from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.builders import maf2poagraph,\
                                             maf2dagmaf,\
                                             dagmaf2poagraph,\
                                             po2poagraph
from pangtreebuild.pangenome.parameters import missings
from pangtreebuild.pangenome import DAGMaf


def build_from_maf(maf: msa.Maf,
                   metadata: Optional[msa.MetadataCSV] = None,
                   datatype: Optional[graph.DataType] = graph.DataType.Nucleotides) -> \
        graph.Poagraph:
    """Builds poagraph from MAF file.

    Args:
        maf: Multialignment as MAF file.
        metadata: Metadata of sequences present in MAF.
        datatype: Type of the processed data (nucleotides/proteins).

    Returns:
        Poagraph based on given input data.
    """

    nodes, sequences = maf2poagraph.get_poagraph(maf, metadata)
    p = graph.Poagraph(nodes, sequences)
    if metadata:
        graph.Poagraph.complement_metadata_for_sequences_absent_in_metadata_provided(p, metadata)
    p.datatype = datatype
    return p


def build_from_dagmaf(maf: msa.Maf,
                      fasta_provider: Optional[missings.FastaProvider] = missings.ConstBaseProvider(
                          missings.MissingBase()),
                      metadata: Optional[msa.MetadataCSV] = None,
                      datatype: Optional[graph.DataType] = graph.DataType.Nucleotides) -> \
        Tuple[graph.Poagraph, DAGMaf.DAGMaf]:
    """Converts MAF to DagMaf and builds poagraph from MAF file.

    Args:
        maf: Multialignment as MAF file.
        fasta_provider: Provider of bases missing in DagMaf.
        metadata: Metadata of sequences present in MAF.
        datatype: Type of the processed data (nucleotides/proteins).

    Returns:
        Tuple: poagraph based on given input data and dagmaf created
            from input MAF.
    """

    dagmaf = maf2dagmaf.get_dagmaf(maf)
    nodes, sequences = dagmaf2poagraph.get_poagraph(dagmaf,
                                                    fasta_provider,
                                                    metadata)
    p = graph.Poagraph(nodes, sequences)
    if metadata:
        graph.Poagraph.complement_metadata_for_sequences_absent_in_metadata_provided(p, metadata)
    p.datatype = datatype
    return p, dagmaf


def build_from_po(po: msa.Po,
                  metadata: Optional[msa.MetadataCSV] = None,
                  datatype: Optional[graph.DataType] = graph.DataType.Nucleotides) -> graph.Poagraph:
    """Builds poagraph from PO file.

    Args:
        po: Multialignment as PO file.
        metadata: Metadata of sequences present in MAF.
        datatype: Type of the processed data (nucleotides/proteins).

    Returns:
        Poagraph based on given input data.
    """

    nodes, sequences = po2poagraph.get_poagraph(po, metadata)
    p = graph.Poagraph(nodes, sequences)
    if metadata:
        graph.Poagraph.complement_metadata_for_sequences_absent_in_metadata_provided(p, metadata)
    p.datatype = datatype
    return p
