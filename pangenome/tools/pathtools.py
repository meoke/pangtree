from pathlib import Path
from io import StringIO


def get_file_content(path: Path) -> str:
    """Returns file content."""

    with open(path) as input_file:
        return input_file.read()


def get_file_content_stringio(path: Path) -> StringIO:
    """Returns file content."""

    return StringIO(get_file_content(path))


def get_child_path(directory: Path, child_name: str) -> Path:
    return directory.joinpath(child_name)


def dir_exists(dir_path):
    return dir_path.exists() and dir_path.is_dir()


def create_dir(cache_dir: Path):
    cache_dir.mkdir()