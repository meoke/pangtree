from data.Poagraph import Poagraph
from tools.cli import get_parser

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    poagraph = Poagraph.build_from_maf(args.multialignment)
    print("h")