import logging

import time

from .userio import cmdargs, pathtools
from .Pangraph import Pangraph


def run_pang(args):
    """Creates Pangraph and runs required algorithms."""
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
    if args.vis:
        p.generate_visualization(pathtools.create_child_dir(args.output, 'vis'))


args = cmdargs.get_validated_args()
logging.info(f'Input arguments: {args}')
start = time.clock()

run_pang(args)

end = time.clock()
processing_time = time.strftime('%H:%M:%S', time.gmtime(end - start))
logging.info(f'DONE! Running time: {processing_time}')