from io import StringIO
from pathlib import Path
import shutil

from datetime import datetime


def create_child_dir(parent_path: Path, child_dir_name: str, add_timestamp: bool = False) -> Path:
    """Creates child dir with given name under parent path."""

    if add_timestamp:
        child_dir_name += _get_current_time()

    child_dir_path = parent_path.joinpath(child_dir_name)
    if not child_dir_path.is_dir():
        child_dir_path.mkdir()
    return child_dir_path.resolve()


def create_default_output_dir(parrent_path: Path) -> Path:
    """Creates timestamped child dir under parent path"""

    output_dir_prefix = 'pang_output'
    current_time = _get_current_time()
    output_dir_name = "_".join([output_dir_prefix, current_time])
    return create_child_dir(parrent_path, output_dir_name)


def remove_dir(dir_path: Path) -> None:
    if dir_path.exists() and dir_path.is_dir():
        shutil.rmtree(dir_path)


def remove_dir_if_empty(dir_path: Path) -> bool:
    """Checks if directory exists and is empty, if yes - removes it."""

    if dir_path.exists() and dir_path.is_dir():
        for item in dir_path.iterdir():
            if not item.is_dir() or (item.is_dir() and [*item.iterdir()]):
                return False
            else:
                item.rmdir()
        dir_path.rmdir()
        return True


def get_child_file_path(directory: Path, file_name: str) -> Path:
    return directory.joinpath(file_name)


def _get_current_time() -> str:
    return datetime.now().strftime('%m_%d__%H_%M_%S')


def get_file_content(path: Path) -> str:
    """Returns file content."""

    with open(path) as input_file:
        return input_file.read()


