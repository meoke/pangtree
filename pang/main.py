import logging.config
import time

from userio import cmdargs, pathtools
from Pangenome import Pangenome
from fileformat.json import writer as jsonwriter


def run_pang(args):
    """Creates Pangraph and runs required algorithms."""

    p = Pangenome(args.multialignment, args.data, as_string=False)
    if args.fasta:
        p.generate_fasta_files(pathtools.create_child_dir(args.output, 'fasta'))
    if args.consensus:
        p.generate_consensus(pathtools.create_child_dir(args.output, 'consensus'),
                             args.consensus,
                             args.hbmin,
                             args.r,
                             args.multiplier,
                             args.stop,
                             args.re_consensus,
                             args.anti_granular
                             )
    if args.vis:
        p.generate_visualization(pathtools.create_child_dir(args.output, 'vis'))
    data_path = pathtools.create_child_dir(args.output, 'data')
    jsonwriter.save(data_path, p)


def cleanup(args: cmdargs.ArgsList)-> None:
    """Removes output directory if it is empty."""

    no_output = pathtools.remove_dir_if_empty(args.output)
    if no_output:
        logging.info("No output was produced.")
    else:
        logging.info(f"Program output is placed in {args.output}")


def main():
    # logging.config.fileConfig('logging.conf')

    args = cmdargs.get_validated_args()
    logging.info(f'Input arguments: {args}')
    start = time.clock()

    try:
        run_pang(args)
    except Exception as e:
        logging.error("Something went wrong...")
        raise e
    finally:
        cleanup(args)

    end = time.clock()
    processing_time = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print(f'DONE! Running time: {processing_time}')


if __name__ == "__main__":
    main()
