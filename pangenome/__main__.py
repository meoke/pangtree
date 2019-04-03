import os
import sys


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangenome')))
from datamodel.input_types import Maf
print(sys.path)
from datamodel.Poagraph import Poagraph
from tools.cli import get_parser
from tools.dependencies import get_fasta_provider


def main():
    parser = get_parser()
    args = parser.parse_args()
    if isinstance(args.multialignment, Maf) and args.raw_maf:
        poagraph = Poagraph.build_from_maf(args.multialignment, args.metadata)
    elif isinstance(args.multialignment, Maf) and not args.raw_maf:
        fasta_provider = get_fasta_provider(args)
        poagraph = Poagraph.build_from_dagmaf(args.multialignment, fasta_provider, args.metadata)
    print(poagraph)


if __name__ == "__main__":
    main()
