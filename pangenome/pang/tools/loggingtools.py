import logging
import logging.config
import os
import sys
from pathlib import Path

from . import pathtools
logging_config_path = Path(os.path.abspath(__file__)).joinpath('../../logging.conf').resolve()
logging.config.fileConfig(logging_config_path)


def add_fileHandler_to_logger(outputdir,
                              logger_name,
                              filename,
                              format="%(levelname)s - %(message)s",
                              propagate=True):
    logger = logging.getLogger(logger_name)
    fh = logging.FileHandler(pathtools.get_child_path(outputdir, filename))
    ft = logging.Formatter(format, datefmt="%x-%X")
    fh.setFormatter(ft)
    logger.propagate = propagate
    logger.addHandler(fh)


def add_consoleHandler_to_logger(logger_name: str, propagate: bool):
    console_logger = logging.getLogger(logger_name)
    console_logger.propagate = propagate
    console_handler = logging.StreamHandler(sys.stdout)
    console_logger.addHandler(console_handler)


def get_logger(loggername):
    return logging.getLogger(loggername)


def get_global_logger():
    return logging.getLogger("")


def disable_all_loggers():
    for logger in [logging.getLogger(name) for name in logging.root.manager.loggerDict]:
        logger.disabled = True

def remove_consoleHandler_from_logger(logger_name):
    logger = logging.getLogger(logger_name)
    for lh in logger.handlers[:]:
        if type(lh) == logging.StreamHandler and lh.stream.name == "<stdout>":
            logger.removeHandler(lh)
