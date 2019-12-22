# Build pangenome and generate affinity trees with different P values.


from pathlib import Path
import os
import sys
from typing import Tuple


def get_relative_path(path_suffix: str) -> Path:
    return Path(os.path.abspath(os.path.join(os.path.dirname(__file__),
                path_suffix)))


from pangtreebuild.pangenome import builder
from pangtreebuild.pangenome.parameters import missings, msa
from pangtreebuild.affinity_tree import parameters as at_params
from pangtreebuild.affinity_tree import builders as at_builders
from pangtreebuild.tools import pathtools
from pangtreebuild.serialization import json, po


p_values = [1, 2]


def run_pangtree(maf_path: Path,
                 fasta_path: Path,
                 output_dir: Path,
                 po_output: bool) -> None:
    output_dir = pathtools.get_child_dir(output_dir,
                                         pathtools.get_current_time())
    print(f"Runing pangtree for maf: {maf_path} and fasta: {fasta_path} "
          f"Output in: {output_dir}, include po file: {po_output}.")

    fasta_provider = missings.FromFile(fasta_path)
    maf = msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path)
    poagraph, dagmaf = builder.build_from_dagmaf(maf, fasta_provider)
    for p in p_values:
        current_output_dir = pathtools.get_child_dir(output_dir, str(p).replace(".", "_"))
        stop = at_params.Stop(0.99)
        at = at_builders.build_affinity_tree(poagraph,
                                             None,
                                             current_output_dir,
                                             stop,
                                             at_params.P(p),
                                             True)

        at_newick = at.as_newick(None, separate_leaves=True)

        pathtools.save_to_file(at_newick,
                               pathtools.get_child_path(current_output_dir, "affinity_tree.newick"))

        if po_output:
            pangenome_po = po.poagraph_to_PangenomePO(poagraph)
            pathtools.save_to_file(pangenome_po,
                                   pathtools.get_child_path(current_output_dir, "poagraph.po"))

        task_params = json.TaskParameters(multialignment_file_path=str(maf_path),
                                          multialignment_format="maf",
                                          datatype="nucleotides",
                                          blosum_file_path="",
                                          output_path=current_output_dir,
                                          fasta_provider=fasta_provider,
                                          fasta_source_file=fasta_path,
                                          consensus_type="tree",
                                          stop=str(stop),
                                          p=str(p),
                                          output_with_nodes=True)
        pangenomejson = json.to_PangenomeJSON(task_parameters=task_params,
                                              poagraph=poagraph,
                                              dagmaf=dagmaf,
                                              affinity_tree=at)

        pangenome_json_str = json.to_json(pangenomejson)
        pathtools.save_to_file(pangenome_json_str,
                               pathtools.get_child_path(current_output_dir, "pangenome.json"))


def read_cmd_args() -> Tuple[Path, Path, Path]:
    if len(sys.argv) != 4:
        print("Execution failed: wrong number of arguments provided. "
              "Provide path to MAF file, P value and "
              "path to the output directory.")
        sys.exit(1)

    maf_path = Path(sys.argv[1])
    fasta_path = Path(sys.argv[2])
    output_path = Path(sys.argv[3])
    return maf_path, fasta_path, output_path


if __name__ == "__main__":
    cmd_line_args = sys.argv
    if cmd_line_args:
        maf_path, fasta_path, output_path = read_cmd_args()
    else:
        maf_path = get_relative_path("../example_data/Simulated/toy_example/f.maf")
        fasta_path = get_relative_path("../example_data/Simulated/toy_example/sequence.fasta")
        output_path = get_relative_path("../output")

    po_output = True
    run_pangtree(maf_path, fasta_path, output_path, po_output)
