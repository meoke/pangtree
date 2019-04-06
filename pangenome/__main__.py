import os
import sys

import consensuses
from consensuses.input_types import Range

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangenome')))
from datamodel.input_types import Maf
print(sys.path)
from datamodel.Poagraph import Poagraph
from tools.cli import get_parser
from tools.dependencies import resolve_fasta_provider
from output.PangenomeJSON import to_PangenomeJSON, TaskParameters, to_json, to_pickle, load_pickle


def main():
    parser = get_parser()
    args = parser.parse_args()
    if isinstance(args.multialignment, Maf) and args.raw_maf:
        poagraph = Poagraph.build_from_maf(args.multialignment, args.metadata)
    elif isinstance(args.multialignment, Maf) and not args.raw_maf:
        fasta_provider = resolve_fasta_provider(args)
        poagraph, dagmaf = Poagraph.build_from_dagmaf(args.multialignment, fasta_provider, args.metadata)
    pangenomejson = to_PangenomeJSON(TaskParameters(), poagraph, dagmaf)
    pagenome_json_str = to_json(pangenomejson)
    print(len(pagenome_json_str))
    pagenome_json_str = to_pickle(pangenomejson)
    print(len(pagenome_json_str))
    pangenom = load_pickle(pagenome_json_str)

    if args.consensus == 'poa':
        consensuses_tree = consensuses.generators.get_poa_consensuses(args.hbmin, Range(args.r))
    elif args.consensus == 'tree':
        consensuses_tree = consensuses.generators.get_consensuses_tree()



if __name__ == "__main__":
    main()
