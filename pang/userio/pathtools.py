from pathlib import Path

from datetime import datetime


def create_child_dir(parent_path: Path, child_dir_name: str) -> Path:
    """Creates child dir with given name under parent path."""

    child_dir_path = parent_path.joinpath(child_dir_name)
    child_dir_path.mkdir()
    return child_dir_path.resolve()


def create_default_output_dir(parrent_path: Path) -> Path:
    """Creates timestamped child dir under parent path"""

    output_dir_prefix = 'pang_output'
    current_time = datetime.now().strftime('%m_%d__%H_%M_%S')
    output_dir_name = "_".join([output_dir_prefix, current_time])
    return create_child_dir(parrent_path, output_dir_name)
