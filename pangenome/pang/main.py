import time
from arguments import cmd_arguments
from tools import pathtools, loggingtools
from Pangenome import Pangenome
from fileformats.json import writer as pangenome_to_json_writer
from arguments.PangenomeParameters import PangenomeParameters

console_logger = loggingtools.get_logger('')


def cleanup(params: PangenomeParameters)-> None:
    """Removes output directory if it is empty."""

    no_output = pathtools.remove_dir_if_empty(params.output_path)
    if no_output:
        console_logger.warning("No output was produced.")
    else:
        console_logger.info(f"Program output is placed in {params.output_path}")


def main():
    program_parameters = cmd_arguments.create_pangenome_parameters()
    if program_parameters.quiet:
        loggingtools.remove_console_handler_from_logger("")

    start = time.clock()
    try:
        pangenome = Pangenome(program_parameters)
        pangenome.run()

        data_path = pathtools.create_child_dir(program_parameters.output_path, 'data')
        pangenome_to_json_writer.save_to_file(data_path, pangenome)
    except Exception as e:
        console_logger.exception(e)
    finally:
        cleanup(program_parameters)

    end = time.clock()
    processing_time = time.strftime('%H:%M:%S', time.gmtime(end - start))
    console_logger.info(f'Running time: {processing_time}')


if __name__ == "__main__":
    main()
