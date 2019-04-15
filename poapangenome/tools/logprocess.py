import logging
import logging.config
import os
import sys
from pathlib import Path

from poapangenome.tools import pathtools

logging_config_path = Path(os.path.abspath(__file__)).joinpath('../../logging.conf').resolve()
logging.config.fileConfig(logging_config_path)


def get_global_logger():
    return logging.getLogger("")


def get_logger(logger_name: str):
    return logging.getLogger(logger_name)


def add_file_handler_to_logger(outputdir,
                               logger_name,
                               filename,
                               handlerformat="%(levelname)s - %(message)s",
                               propagate=True):
    logger = logging.getLogger(logger_name)
    fh = logging.FileHandler(pathtools.get_child_path(outputdir, filename))
    ft = logging.Formatter(handlerformat, datefmt="%x-%X")
    fh.setFormatter(ft)
    logger.propagate = propagate
    logger.addHandler(fh)


def add_console_handler_to_logger(logger_name: str, propagate: bool):
    console_logger = logging.getLogger(logger_name)
    console_logger.propagate = propagate
    console_handler = logging.StreamHandler(sys.stdout)
    console_logger.addHandler(console_handler)


def disable_all_loggers():
    logging.root.disabled = True


def remove_console_handler_from_logger(logger_name):
    logger = logging.getLogger(logger_name)
    for lh in logger.handlers[:]:
        if type(lh) == logging.StreamHandler and lh.stream.name == "<stdout>":
            logger.removeHandler(lh)
