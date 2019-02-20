import logging.config
import time
from userio import cmdargs, pathtools
from Pangenome import Pangenome
from fileformats.json import writer as pangenome_to_json_writer
from userio.PangenomeParameters import PangenomeParameters


def cleanup(params: PangenomeParameters)-> None:
    """Removes output directory if it is empty."""

    no_output = pathtools.remove_dir_if_empty(params.output_path)
    if no_output:
        logging.info("No output was produced.")
    else:
        logging.info(f"Program output is placed in {params.output_path}")


def main():
    program_parameters = cmdargs.create_pangenome_parameters()
    logging.info(f'Program arguments: {str(program_parameters)}')
    start = time.clock()

    try:
        pangenome = Pangenome(program_parameters)
        pangenome.run()
        # todo: result = pangenome.run(program_parameters)

        data_path = pathtools.create_child_dir(program_parameters.output_path, 'data')
        pangenome_to_json_writer.save_to_file(data_path, pangenome)
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
