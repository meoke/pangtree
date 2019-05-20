import os
import shutil
from datetime import datetime
from pathlib import Path
from io import StringIO
from typing import Optional


def get_file_content(path: Path) -> str:
    """Returns file content."""

    with open(path) as input_file:
        return input_file.read()


def get_file_content_stringio(path: Path) -> StringIO:
    """Returns file content."""

    return StringIO(get_file_content(path))


def get_child_path(directory: Path, child_name: str) -> Path:
    return directory.joinpath(child_name)


def dir_exists(dir_path: Path):
    return dir_path.exists() and dir_path.is_dir()


def file_exists(file_path: Path):
    return file_path.is_file()


def create_dir(dir_path: Path):
    dir_path.mkdir()


def get_cwd() -> Path:
    """Returns current working directory."""

    return Path(os.getcwd())


def get_current_time() -> str:
    """Returns current date and time in format MM_DD__HH_MM_SS"""

    return datetime.now().strftime('%m_%d__%H_%M_%S')


def save_to_file(filecontent: str, filename: Path, mode: Optional[str] = 'w') -> None:
    """Saves string to file."""

    with open(filename, mode) as output:
        output.write(filecontent)


def get_child_dir(parent_path: Path, child_dir_name: str):
    child_dir_path = parent_path.joinpath(child_dir_name)
    if not child_dir_path.is_dir():
        child_dir_path.mkdir()
    return child_dir_path.resolve()


def remove_dir(dir_path: Path) -> None:
    if dir_path.exists() and dir_path.is_dir():
        shutil.rmtree(dir_path)
