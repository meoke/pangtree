"""This script is prepared for detaild Ebola virus analysis. Steps:
- build Poagraph based on Ebola multialignment (data/Ebola/genome_whole/input/multialignment.maf)
- run Consensus Tree algorithm (P=0.25, STOP=0.99)
- for specific groups of generated consensuses (Group 1: 1, 2, 3; Group 2: 4, 5, 6, 7, 8; Group 3: sons of 4) calculate compatibility for each group in this group

"""

import os, sys
from pathlib import Path
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
import numpy as np

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
from pangtreebuild.consensus.ConsensusTree import ConsensusTree


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


def get_sim1_consensus_tree(p: float, stop: float, output_dir_name: str) -> Tuple[Poagraph, ConsensusTree]:
    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, output_dir_name)
    consensus_output_dir = pathtools.get_child_dir(output_dir_path, "consensus")
    multialignment_path = current_path.parent.joinpath("../data/Simulated1/input/multialignment.maf")
    metadata_path = current_path.parent.joinpath("../data/Simulated1/input/metadata.csv")
    fasta_path = current_path.parent.joinpath("../data/Simulated1/input/data.fasta")
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
                        fasta_provider='FromFile',
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

    fasta_provider = fp_file.FromFile(fastas_file=fasta_path)

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


def global_compatibilities_analysis(consensus_tree: ConsensusTree, groups: List[List[int]]) -> None:
    def calc_comp(consensuses_group: List[int]):
        pairwise_compatibilities = {}
        for i in consensuses_group:
            i_nodes = set(consensus_tree.nodes[i].consensus_path)
            for j in consensuses_group:
                j_nodes = set(consensus_tree.nodes[j].consensus_path)
                pairwise_compatibilities[f"comp_{i}_{j}"] = len(i_nodes.intersection(j_nodes)) / len(i_nodes)

        for k, v in pairwise_compatibilities.items():
            print(k, v)

    for g in groups:
        calc_comp(g)


def local_compatibilities_analysis(poagraph: Poagraph, consensus_tree: ConsensusTree, groups: List[List[int]]) -> None:
    def produce_chart(x, ys, labels, chart_path):
        fig, ax = plt.subplots()
        for i, y in enumerate(ys):
            ax.plot(x, y, label=labels[i] )

        ax.set(xlabel='POA graph columns IDs', ylabel='Local compatibility to other consensus', title=f"Base consensus: {labels[-1]}")
        ax.legend(loc=1)
        ax.grid()

        fig.savefig(chart_path)

    def produce_local_compatibility_chart(consensuses_group: List[int], column_to_nodes: Dict[ColumnID, NodeID]):
        columns_count = max(column_to_nodes.keys())
        frame_size = 1000
        frame_step = 500
        frame_start = 0

        x = list(range(frame_start, columns_count, frame_step))

        for consensus in consensuses_group:
            chart_path = output_dir_path.joinpath(f"{consensus}.png")
            ys = []
            labels = []
            for consensus_to_compare in set(consensuses_group) - {consensus}:
                y = []

                frame_start = 0
                frame_end = frame_start + frame_size

                while frame_start <= columns_count:
                    frame_columns_ids = range(frame_start, frame_end)
                    frame_nodes = set([node_id for col_id in frame_columns_ids for node_id in column_to_nodes[col_id]])
                    c_nodes_in_frame = frame_nodes.intersection(set(consensus_tree.nodes[consensus].consensus_path))
                    consensus_to_compare_nodes_in_frame = frame_nodes.intersection(set(consensus_tree.nodes[consensus_to_compare].consensus_path))
                    if len(c_nodes_in_frame) != 0:
                        comp = len(c_nodes_in_frame.intersection(consensus_to_compare_nodes_in_frame)) / len(c_nodes_in_frame)
                    else:
                        if len(consensus_to_compare_nodes_in_frame) == 0:
                            comp = 1
                        else:
                            comp = 0
                    y.append(comp)
                    frame_start += frame_step
                    frame_end = min(frame_end + frame_step, columns_count)

                ys.append(y)
                labels.append(str(consensus_to_compare))
            labels.append(consensus)
            produce_chart(x, ys, labels, chart_path)

    column_to_nodes = {node.column_id: [] for node in poagraph.nodes}
    for node in poagraph.nodes:
        column_to_nodes[node.column_id].append(node.node_id)

    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, "charts")

    for g in groups:
        produce_local_compatibility_chart(g, column_to_nodes)

ebola_a = [1, 2, 3]
ebola_b = [4, 5, 6, 7, 8]
ebola_c = [49,50]

sim = [1, 2, 6]

ebola_poagraph, ebola_consensus_tree = get_ebola_consensus_tree(p=0.25, stop=0.99, output_dir_name="output_ebola")
# sim1_poagraph, sim1_consensus_tree = get_sim1_consensus_tree(p=1, stop=0.99, output_dir_name="output_sim1")
# global_compatibilities_analysis(ebola_consensus_tree, [a, b, c])
# local_compatibilities_analysis(sim1_poagraph, sim1_consensus_tree, [sim])
local_compatibilities_analysis(ebola_poagraph, ebola_consensus_tree, [ebola_a, ebola_b, ebola_c])



