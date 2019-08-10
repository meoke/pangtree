"""This script is prepared for detaild Ebola virus analysis - and for preparing Ebola tree of life. Steps:
- build Poagraph based on Ebola multialignment (data/Ebola/genome_whole/input/multialignment.maf)
- run Consensus Tree algorithm (P=0.25, STOP=0.99)
- for specific groups of generated consensuses (Group 1: 1, 2, 3; Group 2: 4, 5, 6, 7, 8; Group 3: sons of 4) calculate compatibility for each group in this group

"""

import os, sys
from pathlib import Path
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))
import pangtreebuild.tools.pathtools as pathtools
import pangtreebuild.datamodel.fasta_providers.FromNCBI as fp_ncbi
import pangtreebuild.datamodel.fasta_providers.FromFile as fp_file
import pangtreebuild.datamodel.input_types as inp
from pangtreebuild.datamodel.Poagraph import Poagraph
from pangtreebuild.datamodel.Node import NodeID, ColumnID
import pangtreebuild.consensus.tree_generator as tree_generator
import pangtreebuild.consensus.input_types as cinp
from pangtreebuild.output.PangenomeJSON import TaskParameters
from pangtreebuild.consensus.cutoffs import MAX2, NODE3
from pangtreebuild.consensus.ConsensusTree import ConsensusTree, ConsensusNode, CompatibilityToPath


def get_ebola_consensus_tree(p: float, stop: float, output_dir_name: str) -> Tuple[Poagraph, ConsensusTree]:
    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, output_dir_name)
    consensus_output_dir = pathtools.get_child_dir(output_dir_path, "consensus")
    multialignment_path = current_path.parent.joinpath("../data/Ebola/genome_whole/input/multialignment.maf")
    metadata_path = current_path.parent.joinpath("../data/Ebola/genome_whole/input/metadata.csv")
    blosum_path = current_path.parent.joinpath("../bin/blosum80.mat")

    tp = TaskParameters(running_time="",
                        multialignment_file_path=multialignment_path,
                        multialignment_format="MAF",
                        datatype="N",
                        metadata_file_path=metadata_path,
                        blosum_file_path=blosum_path,
                        output_path=output_dir_path,
                        output_po=False,
                        output_fasta=False,
                        output_with_nodes=False,
                        verbose=False,
                        raw_maf=False,
                        fasta_provider='FromNCBI',
                        cache=True,
                        missing_base_symbol="",
                        fasta_source_file=None,
                        consensus_type="",
                        hbmin=0.8,
                        max_cutoff_option="MAX2",
                        search_range=None,
                        node_cutoff_option="NODE3",
                        multiplier=None,
                        stop=stop,
                        p=p)

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


def add_leaves(consensus_tree:ConsensusTree):
    new_nodes = []
    for node in consensus_tree.nodes:
        if len(node.children_nodes_ids) == 0 and len(node.sequences_ids) > 1:
            for seq_id in node.sequences_ids:
                consensus_node_id = len(consensus_tree.nodes)+len(new_nodes)
                node.children_nodes_ids.append(consensus_node_id)
                new_nodes.append(ConsensusNode(consensus_id=consensus_node_id,
                                           parent_node_id=node.consensus_id,
                                           children_nodes_ids=[],
                                           sequences_ids=[seq_id],
                                               mincomp=CompatibilityToPath(1.0)

                                           ))

    consensus_tree.nodes.extend(new_nodes)
    return consensus_tree

ebola_poagraph, ebola_consensus_tree = get_ebola_consensus_tree(p=0.25, stop=0.99, output_dir_name="output_ebola")
seq_id_to_metadata = {seq_id: seq.seqmetadata for seq_id, seq in ebola_poagraph.sequences.items()}
ebola_consensus_tree = add_leaves(ebola_consensus_tree)
newick_consensus_tree = ebola_consensus_tree.as_newick(seq_id_to_metadata)

pathtools.save_to_file(newick_consensus_tree, pathtools.get_child_path(Path("ebola_tree_of_life"), "consensus_tree.newick"))


