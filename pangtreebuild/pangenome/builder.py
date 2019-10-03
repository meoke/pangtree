from typing import Optional, Tuple

from pangtreebuild.pangenome.parameters import multialignment
from pangtreebuild.pangenome import poagraph
from pangtreebuild.pangenome.builders import maf2poagraph, maf2dagmaf, dagmaf2poagraph, po2poagraph
from pangtreebuild.pangenome.parameters import missings
from pangtreebuild.pangenome.DAGMaf import DAGMaf #TODO change dagmaf improting


def build_from_maf(maf: multialignment.Maf,
                   metadata: Optional[multialignment.MetadataCSV] = None,
                   datatype: Optional[poagraph.DataType] = poagraph.DataType.Nucleotides) -> poagraph.Poagraph:
    nodes, sequences = maf2poagraph.get_poagraph(maf, metadata)
    p = poagraph.Poagraph(nodes, sequences)
    if metadata:
        poagraph.Poagraph._complement_metadata_for_sequences_absent_in_metadata_provided(p, metadata)
    p.datatype = datatype
    return p


def build_from_dagmaf(maf: multialignment.Maf,
                      fasta_provider: Optional[missings.FastaProvider] = missings.ConstBaseProvider(
                          missings.MissingBase()),
                      metadata: Optional[multialignment.MetadataCSV] = None,
                      datatype: Optional[poagraph.DataType] = poagraph.DataType.Nucleotides) -> Tuple[poagraph.Poagraph, DAGMaf]:
    dagmaf = maf2dagmaf.get_dagmaf(maf)
    nodes, sequences = dagmaf2poagraph.get_poagraph(dagmaf, fasta_provider, metadata)
    p = poagraph.Poagraph(nodes, sequences)
    if metadata:
        poagraph.Poagraph._complement_metadata_for_sequences_absent_in_metadata_provided(p, metadata)
    p.datatype = datatype
    return p, dagmaf


def build_from_po(po: multialignment.Po,
                  metadata: Optional[multialignment.MetadataCSV] = None,
                  datatype: Optional[poagraph.DataType] = poagraph.DataType.Nucleotides) -> poagraph.Poagraph:
    nodes, sequences = po2poagraph.get_poagraph(po, metadata)
    p = poagraph.Poagraph(nodes, sequences)
    if metadata:
        poagraph.Poagraph._complement_metadata_for_sequences_absent_in_metadata_provided(p, metadata)
    p.datatype = datatype
    return p