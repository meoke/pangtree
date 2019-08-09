"""This script is prepared for detaild Ebola virus analysis. Steps:
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


def local_compatibilities_analysis_poagraph_coordinates(poagraph: Poagraph, consensus_tree: ConsensusTree, groups: List[List[int]]) -> None:
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


def local_compatibilities_analysis_consensus_coordinates(poagraph: Poagraph, consensus_tree: ConsensusTree, groups: List[List[int]]) -> None:
    def produce_chart(x, ys, labels, chart_path):
        fig, ax = plt.subplots()
        for i, y in enumerate(ys):
            ax.plot(x, y, label=labels[i] )
        for r in [[469, 2689], [3128,4151], [4478,5459], [6038,8069], [6038,7133], [6038,6933], [8508, 9375], [10344, 11100], [11580,18219]]:
            ax.plot([r[0], r[1]], [1,1  ], color="green")
        ax.set(xlabel='POA graph columns IDs', ylabel='Local compatibility to other consensus', title=f"Base consensus: {labels[-1]}")
        ax.legend(loc=4)
        ax.grid()

        fig.savefig(chart_path)

    class Chart:
        def __init__(self, x, consensus, ys, labels):
            self.x = x
            self.consensus = consensus
            self.ys = ys
            self.labels = labels

    def produce_joint_chart(chart_datas: List[Chart], chart_path):
        fig, axs = plt.subplots(len(chart_datas)+2, 1)
        if len(chart_datas) == 5:
            fig.set_size_inches(14.5, 8)
        elif len(chart_datas) == 3:
            fig.set_size_inches(18.5, 6.5)
        line_objects = []
        line_labels = []
        for i, cd in enumerate(chart_datas):
            for j, y in enumerate(cd.ys):
                # if j == 0 and i == 0:
                axs[len(chart_datas)].plot(cd.x, y, color = 'white')
                c_label = ebola_consensus_labels[cd.labels[j]][0]
                lo = axs[i].plot(cd.x, y, label=c_label, color = ebola_consensus_labels[cd.labels[j]][1])
                line_objects.append(lo)
                if c_label not in line_labels:

                    line_labels.append(c_label)
            axs[i].set_xlabel(f'{ebola_consensus_labels[str(cd.consensus)][0]}')
            axs[i].set_ylim(0, 1)

        fig.legend(line_objects,  # The line objects
                   labels=line_labels,  # The labels for each line
                   loc="lower center",  # Position of legend
                   borderaxespad=0.1,  # Small spacing around legend box
                   title="Legend Title",  # Title for the legend
                    # bbox_to_anchor=(1.1, 1.05)
                   ncol=2
                   )


        for r in [(469, 2689, 1, "NP"),
                  (3128,4151, 1, "VP35"),
                  (4478,5459, 1, "VP40"),
                  (6038,8069, 1, "GP"),
                  (6038,7133, 2, "ssGP"),
                  (6038,6933, 3, "sGP"),
                  (8508, 9375, 1, "VP30"),
                  (10344, 11100, 1, "VP24"),
                  (11580, 18219, 1, "L")]:
        # for r in [(2, 5, 1, "jeden"), (6,10, 1, "dwa")]:
            axs[len(chart_datas)].plot([r[0], r[1]], [r[2], r[2]], color="green")

            axs[len(chart_datas)].annotate(r[3], (r[0], r[2]+0.1))
            # axs[len(chart_datas)].set_xlim(0, 19000)
            # axs[len(chart_datas)].set_xlim(0, 19000)

        for k in [len(chart_datas), len(chart_datas)+1]:
            axs[k].tick_params(
                axis='x',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom=False,  # ticks along the bottom edge are off
                top=False,  # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off
            axs[k].tick_params(
                axis='y',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                left=False,  # ticks along the bottom edge are off
                right=False,  # ticks along the top edge are off
                labelleft=False)  # labels along the bottom edge are off
            axs[k].spines['top'].set_visible(False)
            axs[k].spines['right'].set_visible(False)
            axs[k].spines['bottom'].set_visible(False)
            axs[k].spines['left'].set_visible(False)


        fig.text(0.005, 0.6, 'Compatibility', ha='center', va='center', rotation='vertical')

        fig.tight_layout()
        fig.savefig(chart_path, dpi=100)


    def produce_local_compatibility_chart(consensuses_group: List[int]):
        frame_size = 200
        # frame_size = 5
        frame_step = 200
        # frame_step = 5
        joint_chart_name = "_".join([str(c) for c in consensuses_group])
        joint_chart_path = output_dir_path.joinpath(f"{joint_chart_name}.png")
        chart_datas = []
        for consensus in consensuses_group:
            consensus_path = consensus_tree.nodes[consensus].consensus_path
            consensus_length = len(consensus_path)
            single_chart_path = output_dir_path.joinpath(f"{consensus}.png")
            ys = []
            labels = []
            for consensus_to_compare in set(consensuses_group) - {consensus}:
                y = []
                consensus_to_compare_path = consensus_tree.nodes[consensus_to_compare].consensus_path
                frame_start = 0
                frame_end = frame_start + frame_size
                x=[]
                while frame_start <= consensus_length:
                    frame_nodes_indexes = range(frame_start, frame_end)
                    frame_nodes = set([consensus_path[node_index] for node_index in frame_nodes_indexes])
                    comp = len(frame_nodes.intersection(consensus_to_compare_path)) / len(frame_nodes)
                    y.append(comp)
                    x.append(frame_start)
                    frame_start += frame_step
                    frame_end = min(frame_end + frame_step, consensus_length)

                ys.append(y)
                labels.append(str(consensus_to_compare))
            # labels.append(str(consensus))
            chart_datas.append(Chart(x, consensus, ys, labels))
        produce_joint_chart(chart_datas, joint_chart_path)

    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, "charts_ebola_200_200")

    for g in groups:
        produce_local_compatibility_chart(g)

ebola_a = [1, 2, 3]
ebola_b = [4, 5, 6, 7, 8]
# ebola_c = [49,50]

ebola_consensus_labels = {"1": ("All but Marburg 1987", "red"),
                          "2": ("Marburg 1987 I", "green"),
                          "3": ("Marburg 1987 II", "blue"),
                          "4": ("Ebola 2014, Zaire (DRC) 1976-7, DRC 2007", "blue"),
                          "5": ("Sudan 1976-9", "goldenrod"),
                          "6": ("Reston 1989-90", "forestgreen"),
                          "7": ("Bundibugyo 2007 I", "lightskyblue"),
                          "8": ("Bundibugyo 2007 II", "blueviolet")
                          }
sim = [1, 2, 6]

ebola_poagraph, ebola_consensus_tree = get_ebola_consensus_tree(p=0.25, stop=0.99, output_dir_name="output_ebola")
# sim1_poagraph, sim1_consensus_tree = get_sim1_consensus_tree(p=1, stop=0.99, output_dir_name="output_sim1")
# global_compatibilities_analysis(ebola_consensus_tree, [a, b, c])
# local_compatibilities_analysis(sim1_poagraph, sim1_consensus_tree, [sim])
local_compatibilities_analysis_consensus_coordinates(ebola_poagraph, ebola_consensus_tree, [ebola_a, ebola_b])
# local_compatibilities_analysis_consensus_coordinates(sim1_poagraph, sim1_consensus_tree, [sim])



