import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangpang')))
from consensus import tree_generator, simple_tree_generator
from datamodel.input_types import Maf
from datamodel.Poagraph import Poagraph
from tools import cli, pathtools, logprocess
from output.PangenomeJSON import to_PangenomeJSON, TaskParameters, to_json, to_pickle, load_pickle
from output.PangenomePO import poagraph_to_PangenomePO


def main():
    parser = cli.get_parser()
    args = parser.parse_args()

    if not args.quiet and args.verbose:
        logprocess.add_file_handler_to_logger(args.output_dir, "details", "details.log", propagate=False)
        logprocess.add_file_handler_to_logger(args.output_dir, "", "details.log", propagate=False)
    if args.quiet:
        logprocess.disable_all_loggers()

    poagraph, dagmaf = None, None
    if isinstance(args.multialignment, Maf) and args.raw_maf:
        poagraph = Poagraph.build_from_maf(args.multialignment, args.metadata)
    elif isinstance(args.multialignment, Maf) and not args.raw_maf:
        fasta_provider = cli.resolve_fasta_provider(args)
        poagraph, dagmaf = Poagraph.build_from_dagmaf(args.multialignment, fasta_provider, args.metadata)

    blosum = args.blosum if args.blosum else cli.get_default_blosum(args.missing_symbol)
    consensus_output_dir = pathtools.get_child_dir(args.output_dir, "consensus")
    consensus_tree = None
    if args.consensus == 'poa':
        consensus_tree = simple_tree_generator.get_simple_consensus_tree(poagraph,
                                                                         blosum,
                                                                         consensus_output_dir,
                                                                         args.hbmin,
                                                                         args.verbose)
    elif args.consensus == 'tree':
        max_strategy = cli.resolve_max_strategy(args)
        node_strategy = cli.resolve_node_strategy(args)
        consensus_tree = tree_generator.get_consensus_tree(poagraph,
                                                           blosum,
                                                           consensus_output_dir,
                                                           args.stop,
                                                           args.p,
                                                           max_strategy,
                                                           node_strategy,
                                                           args.verbose)

    if args.output_po:
        pangenome_po = poagraph_to_PangenomePO(poagraph)
        pathtools.save_to_file(pangenome_po, pathtools.get_child_path(args.output_dir, "poagraph.po"))

    pangenomejson = to_PangenomeJSON(TaskParameters(), poagraph, dagmaf, consensus_tree)
    pagenome_json_str = to_json(pangenomejson)
    pathtools.save_to_file(pagenome_json_str, pathtools.get_child_path(args.output_dir, "pangenome.json"))
    # pagenome_pickle_str = to_pickle(pangenomejson)


if __name__ == "__main__":
    main()
