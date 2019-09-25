"""This script is prepared for detaild Ebola virus analysis - and for preparing Ebola tree of life. Steps:
- build Poagraph based on Ebola multialignment (data/Ebola/genome_whole/input/multialignment.maf)
- run Consensus Tree algorithm (P=0.25, STOP=0.99)
- for specific groups of generated affinitytree (Group 1: 1, 2, 3; Group 2: 4, 5, 6, 7, 8; Group 3: sons of 4) calculate compatibility for each group in this group

"""

import os, sys
from pathlib import Path
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt

from pangtreebuild.affinity_tree import simple_tree_generator

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))
import pangtreebuild.tools.pathtools as pathtools
import pangtreebuild.datamodel.fasta_providers.FromNCBI as fp_ncbi
import pangtreebuild.datamodel.fasta_providers.FromFile as fp_file
import pangtreebuild.datamodel.input_types as inp
from pangtreebuild.datamodel.Poagraph import Poagraph
from pangtreebuild.datamodel.Node import NodeID, ColumnID
import pangtreebuild.affinity_tree.tree_generator as tree_generator
import pangtreebuild.affinity_tree.input_types as cinp
from pangtreebuild.output.PangenomeJSON import TaskParameters
from pangtreebuild.affinity_tree.cutoffs import MAX2, NODE3
from pangtreebuild.affinity_tree.ConsensusTree import AffinityTree, AffinityNode, Compatibility


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

    return poagraph, tree_generator.get_affinity_tree(poagraph,
                                                      blosum,
                                                      consensus_output_dir,
                                                      cinp.Stop(stop),
                                                      cinp.P(p),
                                                      MAX2(),
                                                      NODE3(),
                                                      False)


def add_leaves(consensus_tree:AffinityTree):
    new_nodes = []
    for node in consensus_tree.nodes:
        if len(node.children) == 0 and len(node.sequences) > 1:
            for seq_id in node.sequences:
                consensus_node_id = len(consensus_tree.nodes)+len(new_nodes)
                node.children.append(consensus_node_id)
                new_nodes.append(AffinityNode(id=consensus_node_id,
                                              parent=node.id,
                                              children=[],
                                              sequences=[seq_id],
                                              mincomp=Compatibility(1.0)

                                              ))

    consensus_tree.nodes.extend(new_nodes)
    return consensus_tree

# ebola_poagraph, ebola_consensus_tree = get_ebola_consensus_tree(p=0.25, stop=0.99, output_dir_name="output_ebola")
ebola_poagraph, ebola_consensus_tree = get_ebola_poa_tree(hbmin=0.8, output_dir_name="output_ebola")

seq_id_to_metadata = {seq_id: seq.seqmetadata for seq_id, seq in ebola_poagraph.sequences.items()}
ebola_consensus_tree = add_leaves(ebola_consensus_tree)
newick_consensus_tree = ebola_consensus_tree.as_newick(seq_id_to_metadata)

pathtools.save_to_file(newick_consensus_tree, pathtools.get_child_path(Path("ebola_poa_tree_of_life"), "consensus_tree_poa.newick"))


