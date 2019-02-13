import logging.config
import time
from io import StringIO

from userio import cmdargs, pathtools
from Pangenome import Pangenome
from fileformats.json import writer as jsonwriter
from userio.ProgramParameters import ProgramParameters, ConsensusAlgorithm
from userio.cmdargs import get_fasta_complementation_option
from userio.pathtools import get_file_content_as_stringio


def run_pang(params: ProgramParameters):
    """Creates Pangraph and runs required algorithms."""

    p = Pangenome(pathtools.get_file_content(params.metadata_file_path))

    multialignment = get_file_content_as_stringio(params.multialignment_file_path)

    if params.not_dag:
        p.build_from_maf(multialignment)
    else:
        fasta_complementation_option = params.fasta_complementation
        p.build_from_maf_firstly_converted_to_dag(multialignment, fasta_complementation_option, params.local_fasta_dirpath)

    if params.generate_fasta:
        p.generate_fasta_files(pathtools.create_child_dir(params.output_path, 'fasta'))

    if params.consensus_type != ConsensusAlgorithm.No:
        p.generate_consensus(pathtools.create_child_dir(params.output_path, 'consensus'),
                             params.consensus_type,
                             params.hbmin,
                             params.r,
                             params.multiplier,
                             params.stop,
                             params.re_consensus,
                             params.anti_granular
                             )
    # if args.vis:
    #     p.generate_visualization(pathtools.create_child_dir(args.output, 'vis'))
    data_path = pathtools.create_child_dir(params.output_path, 'data')
    jsonwriter.save(data_path, p, params)


def cleanup(params: ProgramParameters)-> None:
    """Removes output directory if it is empty."""

    no_output = pathtools.remove_dir_if_empty(params.output_path)
    if no_output:
        logging.info("No output was produced.")
    else:
        logging.info(f"Program output is placed in {params.output_path}")


def main():
    # logging.config.fileConfig('logging.conf')

    program_parameters = cmdargs.get_validated_args()
    logging.info(f'Program arguments: {str(program_parameters)}')
    start = time.clock()

    try:
        run_pang(program_parameters)
    except Exception as e:
        logging.error("Something went wrong...")
        raise e
    finally:
        cleanup(program_parameters)

    end = time.clock()
    processing_time = time.strftime('%H:%M:%S', time.gmtime(end - start))
    print(f'DONE! Running time: {processing_time}')


if __name__ == "__main__":
    main()
