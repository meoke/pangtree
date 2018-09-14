from io import cmdargs, pathtools
from pangraph.pangraph import Pangraph

args = cmdargs.get_validated_args()
p = Pangraph(args.multialignment, args.data)

if args.fasta:
    p.generate_fasta_files(pathtools.create_child_dir(args.output, 'fasta'))
if args.consensus:
    p.generate_consensus(pathtools.create_child_dir(args.output, 'consensus'),
                         args.hbmin,
                         args.mincomp,
                         args.r,
                         args.multiplier,
                         args.stop,
                         args.re_consensus)
if args.viz:
    p.generate_visualization(pathtools.create_child_dir(args.output, 'vis'))
#Kocham Cie!