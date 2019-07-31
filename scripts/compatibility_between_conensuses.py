"""This script is prepared for detaild Ebola virus analysis. Steps:
- build Poagraph based on Ebola multialignment (data/Ebola/genome_whole/input/multialignment.maf)
- run Consensus Tree algorithm (P=0.25, STOP=0.99)
- for specific groups of generated consensuses (Group 1: 1, 2, 3; Group 2: 4, 5, 6, 7, 8; Group 3: sons of 4) calculate compatibility for each group in this group

"""

import os, sys
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))
import pangtreebuild.tools.pathtools as pathtools
import pangtreebuild.datamodel.fasta_providers.FromNCBI as fp_ncbi
import pangtreebuild.datamodel.input_types as inp
from pangtreebuild.datamodel.Poagraph import Poagraph
import pangtreebuild.consensus.tree_generator as tree_generator
import pangtreebuild.consensus.input_types as cinp
from pangtreebuild.output.PangenomeJSON import to_PangenomeJSON, to_json, TaskParameters
from pangtreebuild.consensus.cutoffs import MAX2, NODE3

current_path = Path(os.path.abspath(__file__)).resolve()
output_dir = pathtools.get_child_dir(current_path.parent, "output1")
consensus_output_dir = pathtools.get_child_dir(output_dir, "consensus2")

tp = TaskParameters(running_time="",
              multialignment_file_path="",
              multialignment_format="MAF",
              datatype="N",
              metadata_file_path="",
              blosum_file_path="",
              output_path="",
              output_po=False,
              output_fasta=False,
              output_with_nodes=False,
              verbose=False,
              raw_maf=False,
              fasta_provider='File',
              cache=True,
              missing_base_symbol="",
              fasta_source_file="",
              consensus_type="",
              hbmin=0.8,
              max_cutoff_option="MAX2",
              search_range=0,
              node_cutoff_option="NODE3",
              multiplier=0,
              stop=0.99,
              p=0.25)

fasta_provider = fp_ncbi.FromNCBI(use_cache=True)

multialignment_path = "../data/Ebola/genome_whole/input/multialignment.maf"
multialignment_content = pathtools.get_file_content_stringio(current_path.parent.joinpath(multialignment_path))
multialignment = inp.Maf(file_content=multialignment_content, filename=multialignment_path)


metadata_path = "../data/Ebola/genome_whole/input/metadata.csv"
metadata_content = pathtools.get_file_content_stringio(current_path.parent.joinpath(metadata_path))
metadata = inp.MetadataCSV(filecontent=metadata_content, filename=metadata_path)

poagraph, dagmaf = Poagraph.build_from_dagmaf(multialignment, fasta_provider, metadata)
blosum_path = current_path.parent.joinpath("../bin/blosum80.mat")
blosum_content = pathtools.get_file_content_stringio(path=blosum_path)
blosum = cinp.Blosum(blosum_content, blosum_path)

consensus_tree = tree_generator.get_consensus_tree(poagraph,
                                                   blosum,
                                                   consensus_output_dir,
                                                   cinp.Stop(0.99),
                                                   cinp.P(0.25),
                                                   MAX2(),
                                                   NODE3(),
                                                   False)


def calc_comp(s):
    s_comp = {}
    for i in s:
        i_nodes = set(consensus_tree.nodes[i].consensus_path)
        for j in s:
            j_nodes = set(consensus_tree.nodes[j].consensus_path)
            s_comp[f"comp_{i}_{j}"] = len(i_nodes.intersection(j_nodes)) / len(i_nodes)

    for k, v in s_comp.items():
        print(k, v)


a = [1, 2, 3]
b = [4, 5, 6, 7, 8]

calc_comp(a)
calc_comp(b)