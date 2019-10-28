"""This script is prepared for detaild Ebola virus analysis. Steps:
- build Poagraph based on Ebola multialignment (data/Ebola/genome_whole/input/multialignment.maf)
- run Consensus Tree algorithm (P=0.25, STOP=0.99)
- for specific groups of generated affinitytree (Group 1: 1, 2, 3; Group 2: 4, 5, 6, 7, 8; Group 3: sons of 4) calculate compatibility for each group in this group

"""

import os, sys
from pathlib import Path
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import matplotlib

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))
import pangtreebuild.tools.pathtools as pathtools
import pangtreebuild.pangenome.fasta_providers.FromNCBI as fp_ncbi
import pangtreebuild.pangenome.fasta_providers.FromFile as fp_file
import pangtreebuild.pangenome.input_types as inp
from pangtreebuild.pangenome.Poagraph import Poagraph
from pangtreebuild.pangenome.Node import NodeID, ColumnID
import pangtreebuild.affinity_tree.tree_generator as tree_generator
import pangtreebuild.affinity_tree.input_types as cinp
from pangtreebuild.output.PangenomeJSON import TaskParameters
from pangtreebuild.affinity_tree.cutoffs import MAX2, NODE3
from pangtreebuild.affinity_tree.ConsensusTree import AffinityTree


def get_ebola_consensus_tree(p: float, stop: float, output_dir_name: str) -> Tuple[Poagraph, AffinityTree]:
    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, output_dir_name)
    consensus_output_dir = pathtools.get_child_dir(output_dir_path, "consensus")
    multialignment_path = current_path.parent.joinpath("../data/Ebola/multialignment.maf")
    metadata_path = current_path.parent.joinpath("../data/Ebola/metadata.csv")
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

    return poagraph, tree_generator.get_affinity_tree(poagraph,
                                                      blosum,
                                                      consensus_output_dir,
                                                      cinp.Stop(stop),
                                                      cinp.P(p),
                                                      MAX2(),
                                                      NODE3(),
                                                      False)


def global_compatibilities_analysis(consensus_tree: AffinityTree, groups: List[List[int]]) -> None:
    def calc_comp(consensuses_group: List[int]):
        pairwise_compatibilities = {}
        for i in consensuses_group:
            i_nodes = set(consensus_tree.nodes[i].consensus)
            for j in consensuses_group:
                j_nodes = set(consensus_tree.nodes[j].consensus)
                pairwise_compatibilities[f"comp_{i}_{j}"] = len(i_nodes.intersection(j_nodes)) / len(i_nodes)

        for k, v in pairwise_compatibilities.items():
            print(k, v)

    for g in groups:
        calc_comp(g)


def local_compatibilities_analysis_consensus_coordinates(poagraph: Poagraph, consensus_tree: AffinityTree, groups: List[List[int]]) -> None:
    def produce_chart(x, ys, labels, chart_path):
        fig, ax = plt.subplots()
        for i, y in enumerate(ys):
            ax.plot(x, y, label=labels[i] )
        for r in [[469, 2689], [3128,4151], [4478,5459], [6038,8069], [6038,7133], [6038,6933], [8508, 9375], [10344, 11100], [11580,18219]]:
            ax.plot([r[0], r[1]], [1,1  ], color="green")
        ax.set(xlabel='POA graph columns IDs', ylabel='Local compatibility', title=f"Base consensus: {labels[-1]}")
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
        matplotlib.rc('text', usetex=True)
        fig, axs = plt.subplots(len(chart_datas)+1, 1)
        if len(chart_datas) == 5:
            fig.set_size_inches(14.5, 8)
        elif len(chart_datas) == 3:
            fig.set_size_inches(18.5, 6.5)

        genes = [(469, 2689, "NP"),
                  (3128,4151, "VP35"),
                  (4478,5459, "VP40"),
                  (6038,8069, "GP"),
                  (8508, 9375, "VP30"),
                  (10344, 11100, "VP24"),
                  (11580, 18219, "L")]

        line_objects = []
        line_labels = []
        for i, cd in enumerate(chart_datas):
            for j, y in enumerate(cd.ys):
                # if j == 0 and i == 0:
                axs[len(chart_datas)].plot(cd.x, [0 for _ in y], color = 'white')
                c_label = ebola_consensus_labels[cd.labels[j]][0]
                color = ebola_consensus_labels[cd.labels[j]][1]
                c_label = r"\textit{" + c_label + "}"
                if c_label in line_labels:
                    c_label = '_' + c_label
                line_labels.append(c_label)

                lo = axs[i].plot(cd.x, y, label=c_label, color = color)
                line_objects.append(lo)

            # plot coding areas
            for g in genes:
                axs[i].add_patch(Rectangle((g[0], 0), g[1]-g[0], 1,
                      color= (0.1, 0.2, 0.5, 0.3)))
            x_label = r"\textit{" + f'{ebola_consensus_labels[str(cd.consensus)][0]}' + "}"
            axs[i].set_xlabel(x_label)
            axs[i].set_ylim(0, 1)

        fig.legend(loc="lower right",  # Position of legend
                   borderaxespad=1,  # Small spacing around legend box
                   title="Consensus sequence",  # Title for the legend
                   ncol=2,
                   )

        for r in [(469, 2689,  "NP"),
                  (3128,4151,  "VP35"),
                  (4478,5459,  "VP40"),
                  (6038,8069,  "GP"),
                  # (6038,7133, 2, "ssGP"),
                  # (6038,6933, 3, "sGP"),
                  (8508, 9375,  "VP30"),
                  (10344, 11100,  "VP24"),
                  (11580, 18219,  "L")]:

            axs[len(chart_datas)].plot([r[0], r[1]], [1, 1], color="white")
            axs[len(chart_datas)].annotate(r[2], (r[0], 1))

        for k in [len(chart_datas), len(chart_datas)]:
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


        fig.text(0.005, 0.6, 'Local compatibility', ha='center', va='center', rotation='vertical')

        fig.tight_layout()
        fig.savefig(chart_path, dpi=100)


    def produce_local_compatibility_chart(consensuses_group: List[int]):
        frame_size = 400
        frame_step = 200
        joint_chart_name = "_".join([str(c) for c in consensuses_group])
        joint_chart_path = output_dir_path.joinpath(f"{joint_chart_name}.png")
        chart_datas = []
        for consensus in consensuses_group:
            consensus_path = consensus_tree.nodes[consensus].consensus
            consensus_length = len(consensus_path)
            ys = []
            labels = []
            for consensus_to_compare in set(consensuses_group) - {consensus}:
                y = []
                consensus_to_compare_path = consensus_tree.nodes[consensus_to_compare].consensus
                frame_start = 0
                frame_end = frame_start + frame_size
                x=[]
                while frame_start <= consensus_length:
                    frame_nodes_indexes = range(frame_start, frame_end)
                    frame_nodes = set([consensus_path[node_index] for node_index in frame_nodes_indexes])
                    comp = len(frame_nodes.intersection(consensus_to_compare_path)) / len(frame_nodes)
                    y.append(comp)
                    x.append(frame_start + frame_step/2)
                    frame_start += frame_step
                    frame_end = min(frame_end + frame_step, consensus_length)

                ys.append(y)
                labels.append(str(consensus_to_compare))
            chart_datas.append(Chart(x, consensus, ys, labels))
        produce_joint_chart(chart_datas, joint_chart_path)

    current_path = Path(os.path.abspath(__file__)).resolve()
    output_dir_path = pathtools.get_child_dir(current_path.parent, "charts_ebola_400_200_middle_prostokaty")

    for g in groups:
        produce_local_compatibility_chart(g)

ebola_a = [1, 2, 3]
ebola_b = [4, 5, 6, 7, 8]

ebola_consensus_labels = {"1": ("Zaire ebolavirus", "sienna"),
                          "2": ("Marburgvirus 1980", "green"),
                          "3": ("Marburgvirus 1987", "orange"),
                          "4": ("Zaire (DRC) ebolavirus", "grey"),
                          "5": ("Sudan ebolavirus", "lightgreen"),
                          "6": ("Reston ebolavirus", "darkviolet"),
                          "7": ("Bundibugyo ebolavirus", "hotpink"),
                          "8": ("Tai Forest ebolavirus", "yellow")
                          }
sim = [1, 2, 6]

os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2019/bin/x86_64-linux'


ebola_poagraph, ebola_consensus_tree = get_ebola_consensus_tree(p=0.25, stop=0.99, output_dir_name="output_ebola")
local_compatibilities_analysis_consensus_coordinates(ebola_poagraph, ebola_consensus_tree, [ebola_a, ebola_b])



