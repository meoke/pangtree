from pathlib import Path
from io import StringIO


def get_file_content(path: Path) -> str:
    """Returns file content."""

    with open(path) as input_file:
        return input_file.read()


def get_file_content_stringio(path: Path) -> StringIO:
    """Returns file content."""

    return StringIO(get_file_content(path))
