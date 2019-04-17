import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../poapangenome')))
from poapangenome.consensus import tree_generator, simple_tree_generator
from poapangenome.datamodel.input_types import Maf, Po
from poapangenome.datamodel.Poagraph import Poagraph
from poapangenome.tools import cli, pathtools, logprocess
from poapangenome.output.PangenomeJSON import to_PangenomeJSON, TaskParameters, to_json, to_pickle, load_pickle, str_to_PangenomeJSON
from poapangenome.output.PangenomePO import poagraph_to_PangenomePO
from poapangenome.output.PangenomeFASTA import poagraph_to_fasta, consensuses_tree_to_fasta


def main():
    parser = cli.get_parser()
    args = parser.parse_args()

    start = time.time()
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
    elif isinstance(args.multialignment, Po):
        poagraph = Poagraph.build_from_po(args.multialignment, args.metadata)

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

    if args.output_fasta:
        sequences_fasta = poagraph_to_fasta(poagraph)
        pathtools.save_to_file(sequences_fasta, pathtools.get_child_path(args.output_dir, "sequences.fasta"))
        if consensus_tree:
            consensuses_fasta = consensuses_tree_to_fasta(poagraph, consensus_tree)
            pathtools.save_to_file(consensuses_fasta, pathtools.get_child_path(args.output_dir, "consensuses.fasta"))

    end = time.time()
    pangenomejson = to_PangenomeJSON(task_parameters=cli.get_task_parameters(args, running_time=f"{end-start}s"),
                                     poagraph=poagraph,
                                     dagmaf=dagmaf,
                                     consensuses_tree=consensus_tree)
    pangenome_json_str = to_json(pangenomejson)
    pangenome_json = str_to_PangenomeJSON(pangenome_json_str)
    pathtools.save_to_file(pangenome_json_str, pathtools.get_child_path(args.output_dir, "pangenome.json"))
    # pagenome_pickle_str = to_pickle(pangenomejson)
    # pathtools.save_to_file(pagenome_pickle_str, pathtools.get_child_path(args.output_dir, "pangenome.pickle"), 'wb')
    # jsonpangenome = load_pickle(pagenome_pickle_str)
    print()

if __name__ == "__main__":
    main()
