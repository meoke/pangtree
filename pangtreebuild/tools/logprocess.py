import logging
import logging.config
import os
import sys
from pathlib import Path

from pangtreebuild.tools import pathtools

logging_config_path = Path(os.path.abspath(__file__)).joinpath('../../logging.conf').resolve()
logging.config.fileConfig(logging_config_path)


def get_global_logger() -> logging.Logger:
    """Returns global logger."""

    return logging.getLogger("")


def get_logger(logger_name: str) -> logging.Logger:
    """Returns logger of given logger_name.

    Args:
        logger_name: Name of the logger to return.

    Returns:
        Logger with name logger_name.
    """

    return logging.getLogger(logger_name)


def add_file_handler_to_logger(outputdir: Path,
                               logger_name: str,
                               filename: str,
                               handlerformat: str = "%(levelname)s - %(message)s",
                               propagate: bool = True) -> None:
    """Adds new file handle to given logger.

    Args:
        outputdir: Directory where the file will be stored.
        logger_name: Logger to modify.
        filename: New logger name.
        handlerformat: Log message format.
        propagate: Whether logger should propagate.
    """

    logger = logging.getLogger(logger_name)
    fh = logging.FileHandler(pathtools.get_child_path(outputdir, filename))
    ft = logging.Formatter(handlerformat, datefmt="%x-%X")
    fh.setFormatter(ft)
    logger.propagate = propagate
    logger.addHandler(fh)


def add_console_handler_to_logger(logger_name: str,
                                  propagate: bool):
    """Adds new console handle to given logger.

    Args:
        logger_name: Logger to modify.
        propagate: Whether logger should propagate.
    """

    console_logger = logging.getLogger(logger_name)
    console_logger.propagate = propagate
    console_handler = logging.StreamHandler(sys.stdout)
    console_logger.addHandler(console_handler)


def disable_all_loggers() -> None:
    """Turns off all loggers."""

    logging.root.disabled = True


def remove_console_handler_from_root_logger() -> None:
    """Turn off console logger completely."""

    logger = logging.getLogger("root").parent
    for lh in logger.handlers[:]:
        if type(lh) == logging.StreamHandler and lh.stream.name == "<stdout>":
            logger.removeHandler(lh)
