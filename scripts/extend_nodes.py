"""This script is prepared for detaild Ebola virus analysis - and for preparing Ebola tree of life. Steps:
- build Poagraph based on Ebola multialignment (data/Ebola/genome_whole/input/multialignment.maf)
- run Consensus Tree algorithm (P=0.25, STOP=0.99)
- for specific groups of generated affinitytree (Group 1: 1, 2, 3; Group 2: 4, 5, 6, 7, 8; Group 3: sons of 4) calculate compatibility for each group in this group

"""

import os, sys
from pathlib import Path
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt

from pangtreebuild.affinitytree import simple_tree_generator
from pangtreebuild.datamodel.Sequence import Sequence, SequenceID

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))
import pangtreebuild.tools.pathtools as pathtools
from pangtreebuild.datamodel.Poagraph import Poagraph
from pangtreebuild.datamodel.Node import Node, NodeID, Base
from pangtreebuild.affinitytree.ConsensusTree import AffinityTree, AffinityNode, Compatibility, ConsensusNodeID
from pangtreebuild.output.PangenomeJSON import str_to_PangenomeJSON, PangenomeJSON


def convert_jsonpangenome(path: Path) -> Tuple[Poagraph, AffinityTree]:
    jsonpangenome = str_to_PangenomeJSON(pathtools.get_file_content(path))
    poagraph = Poagraph(nodes=[Node(node_id=NodeID(n.id), base=Base(n.base)) for n in jsonpangenome.nodes],
                    sequences={SequenceID(seq.sequence_str_id): Sequence(seqid=SequenceID(seq.sequence_str_id), paths=[], seqmetadata=seq.metadata) for seq in jsonpangenome.sequences},
                    )
    seq_int_to_str_ids = {s.sequence_int_id: s.sequence_str_id for s in jsonpangenome.sequences}
    consensus_tree = AffinityTree()
    consensus_tree.nodes = [AffinityNode(id=ConsensusNodeID(c.consensus_node_id),
                                         parent=ConsensusNodeID(c.parent),
                                         children=[ConsensusNodeID(child) for child in c.children],
                                         sequences=[SequenceID(seq_int_to_str_ids[seq_int_id]) for seq_int_id in c.sequences_int_ids],
                                         mincomp=Compatibility(c.mincomp)) for c in jsonpangenome.affinitytree]

    return poagraph, consensus_tree

def get_ebola_poa_tree(hbmin: float, output_dir_name: str) -> Tuple[Poagraph, AffinityTree]:
    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, output_dir_name)
    consensus_output_dir = pathtools.get_child_dir(output_dir_path, "consensus")
    multialignment_path = current_path.parent.joinpath("../data/Ebola/genome_whole/input/multialignment.maf")
    metadata_path = current_path.parent.joinpath("../data/Ebola/genome_whole/input/metadata.csv")
    blosum_path = current_path.parent.joinpath("../bin/blosum80.mat")

    fasta_provider = fp_ncbi.FromNCBI(use_cache=True)

    multialignment_content = pathtools.get_file_content_stringio(multialignment_path)
    multialignment = inp.Maf(file_content=multialignment_content, filename=multialignment_path)

    metadata_content = pathtools.get_file_content_stringio(metadata_path)
    metadata = inp.MetadataCSV(filecontent=metadata_content, filename=metadata_path)

    poagraph, dagmaf = Poagraph.build_from_dagmaf(multialignment, fasta_provider, metadata)

    blosum_content = pathtools.get_file_content_stringio(path=blosum_path)
    blosum = cinp.Blosum(blosum_content, blosum_path)

    return poagraph, simple_tree_generator.get_simple_affinity_tree(poagraph,
                                                                    blosum,
                                                                    consensus_output_dir,
                                                                    cinp.Hbmin(hbmin),
                                                                    False)

def get_ebola_consensus_tree(p: float, stop: float, output_dir_name: str) -> Tuple[Poagraph, AffinityTree]:
    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, output_dir_name)
    consensus_output_dir = pathtools.get_child_dir(output_dir_path, "consensus")
    multialignment_path = current_path.parent.joinpath("../data/Ebola/genome_whole/input/multialignment.maf")
    metadata_path = current_path.parent.joinpath("../data/Ebola/genome_whole/input/metadata.csv")
    blosum_path = current_path.parent.joinpath("../bin/blosum80.mat")


    fasta_provider = fp_ncbi.FromNCBI(use_cache=True)

    multialignment_content = pathtools.get_file_content_stringio(multialignment_path)
    multialignment = inp.Maf(file_content=multialignment_content, filename=multialignment_path)

    metadata_content = pathtools.get_file_content_stringio(metadata_path)
    metadata = inp.MetadataCSV(filecontent=metadata_content, filename=metadata_path)

    poagraph, dagmaf = Poagraph.build_from_dagmaf(multialignment, fasta_provider, metadata)

    blosum_content = pathtools.get_file_content_stringio(path=blosum_path)
    blosum = cinp.Blosum(blosum_content, blosum_path)

    return poagraph, tree_generator.get_consensus_tree(poagraph,
                                             blosum,
                                             consensus_output_dir,
                                             cinp.Stop(stop),
                                             cinp.P(p),
                                             MAX2(),
                                             NODE3(),
                                             False)


for sample in range(0, 10):
    for p in ["025", "05", "1", "2", "4"]:
        dir_path = f"../data/Simulated3_fewcycles/output/analysis_1108_040134/{sample}_{p}__output/"
        # dir_path = f"../data/Simulated1/output/analysis_0808_035057/{sample}_{p}__output/"
        poagraph, consensus_tree = convert_jsonpangenome(path=Path(f"{dir_path}/pangenome.json"))
        seq_id_to_metadata = {seq_id: seq.seqmetadata for seq_id, seq in poagraph.sequences.items()}
        seq_id_to_metadata = None
        newick_consensus_tree = consensus_tree.as_newick(seq_id_to_metadata, separate_leaves=True)
        pathtools.save_to_file(newick_consensus_tree, Path(f"{dir_path}/consensus_tree_extended.newick"))


