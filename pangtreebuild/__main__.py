import os
import sys
import datetime


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangtreebuild')))
from pangtreebuild.affinity_tree import builders as at_builders
from pangtreebuild.output import fasta
from pangtreebuild.output import json
from pangtreebuild.output import po
from pangtreebuild.pangenome import builder
from pangtreebuild.pangenome.parameters import missings
from pangtreebuild.pangenome.parameters import multialignment
from pangtreebuild.tools import cli, pathtools, logprocess


def main():
    parser = cli.get_parser()
    args = parser.parse_args()
    start = datetime.datetime.now()
    if not args.quiet and args.verbose:
        logprocess.add_file_handler_to_logger(args.output_dir, "details", "details.log", propagate=False)
        logprocess.add_file_handler_to_logger(args.output_dir, "", "details.log", propagate=False)
    if args.quiet:
        logprocess.disable_all_loggers()

    poagraph, dagmaf, fasta_provider = None, None, None
    if isinstance(args.multialignment, multialignment.Maf) and args.raw_maf:
        poagraph = builder.build_from_maf(args.multialignment, args.metadata)
    elif isinstance(args.multialignment, multialignment.Maf) and not args.raw_maf:
        fasta_provider = cli.resolve_fasta_provider(args)
        poagraph, dagmaf = builder.build_from_dagmaf(args.multialignment, fasta_provider, args.metadata)
    elif isinstance(args.multialignment, multialignment.Po):
        poagraph = builder.build_from_po(args.multialignment, args.metadata)

    affinity_tree = None
    if args.consensus is not None:
        blosum = args.blosum if args.blosum else cli.get_default_blosum()
        if fasta_provider is not None and isinstance(fasta_provider, missings.ConstBaseProvider):
            blosum.check_if_symbol_is_present(fasta_provider.missing_base.as_str())

        consensus_output_dir = pathtools.get_child_dir(args.output_dir, "consensus")

        if args.consensus == 'poa':
            affinity_tree = at_builders.build_poa_affinity_tree(poagraph,
                                                                            blosum,
                                                                            consensus_output_dir,
                                                                            args.hbmin,
                                                                            args.verbose)
        elif args.consensus == 'tree':
            affinity_tree = at_builders.build_affinity_tree(poagraph,
                                                              blosum,
                                                              consensus_output_dir,
                                                              args.stop,
                                                              args.p,
                                                              args.verbose)

        if args.metadata is not None:
            seq_id_to_metadata = {seq_id: seq.seqmetadata for seq_id, seq in poagraph.sequences.items()}
        else:
            seq_id_to_metadata = None

        newick_affinity_tree = affinity_tree.as_newick(seq_id_to_metadata, separate_leaves=False)
        newick_affinity_tree_extended = affinity_tree.as_newick(seq_id_to_metadata, separate_leaves=True)

        pathtools.save_to_file(newick_affinity_tree, pathtools.get_child_path(args.output_dir, "affinity_tree.newick"))
        pathtools.save_to_file(newick_affinity_tree_extended, pathtools.get_child_path(args.output_dir, "affinity_tree_extended.newick"))


    if args.output_po:
        pangenome_po = po.poagraph_to_PangenomePO(poagraph)
        pathtools.save_to_file(pangenome_po, pathtools.get_child_path(args.output_dir, "poagraph.po"))

    if args.output_fasta:
        sequences_fasta = fasta.poagraph_to_fasta(poagraph)
        pathtools.save_to_file(sequences_fasta, pathtools.get_child_path(args.output_dir, "_sequences.fasta"))
        if affinity_tree:
            consensuses_fasta = fasta.affinity_tree_to_fasta(poagraph, affinity_tree)
            pathtools.save_to_file(consensuses_fasta, pathtools.get_child_path(args.output_dir, "affinitytree.fasta"))

    end = datetime.datetime.now()
    pangenomejson = json.to_PangenomeJSON(task_parameters=cli.get_task_parameters(args, running_time=f"{end-start}s"),
                                     poagraph=poagraph,
                                     dagmaf=dagmaf,
                                     affinity_tree=affinity_tree)

    pangenome_json_str = json.to_json(pangenomejson)
    pathtools.save_to_file(pangenome_json_str, pathtools.get_child_path(args.output_dir, "pangenome.json"))


if __name__ == "__main__":
    main()
