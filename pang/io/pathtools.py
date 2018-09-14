from pathlib import Path


def create_child_dir(parent_path, child_dir_name):
    child_dir_path = Path(parent_path).joinpath(child_dir_name)
    child_dir_path.mkdir()
    return child_dir_path.resolve()