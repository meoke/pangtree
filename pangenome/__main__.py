import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../pangenome')))
print(sys.path)
from datamodel.Poagraph import Poagraph
from tools.cli import get_parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    poagraph = Poagraph.build_from_maf(args.multialignment)
    print(poagraph)


if __name__ == "__main__":
    main()
