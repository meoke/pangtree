#import os
from pathlib import Path
import shutil

def get_file_name_without_extension(file_path):
    return Path(file_path).stem


def get_parentdir_name(file_path):
    return Path(file_path).parent


def is_dir(path):
    return Path(path).is_dir()


def create_next_sibling_dir(path, sibling_dir_prefix):
    parent_path = get_parentdir_name(path)
    existing_prefixed_dirs = sorted(list(parent_path.glob("*"+sibling_dir_prefix+"*")))
    if existing_prefixed_dirs:
        try:
            last_dir_suffix = int((existing_prefixed_dirs[-1].stem).split("_")[-1])
        except:
            last_dir_suffix = -1
    else:
        last_dir_suffix = -1
    new_dir_suffix_value = last_dir_suffix + 1
    new_dir_suffix = str(new_dir_suffix_value) if new_dir_suffix_value > 9 else '0' + str(new_dir_suffix_value)
    new_dir_name = sibling_dir_prefix + '_' + new_dir_suffix
    new_dir_path = parent_path.joinpath(new_dir_name)
    new_dir_path.mkdir()
    return new_dir_path.resolve()

def create_child_dir(parent_path, dir_name):
    child_dir_path = Path(parent_path).joinpath(dir_name)
    child_dir_path.mkdir()
    return child_dir_path.resolve()


def copy_file_to_dir(file_path, destination_dir):
    return shutil.copy(str(Path(file_path)), str(Path(destination_dir)))


def remove_dir(path):
    for i in Path(path).iterdir():
        if is_dir(i):
            remove_dir(i)
        else:
            i.unlink()
    Path(path).rmdir()


def save_text(text, path, file_name):
    file_path = Path(path).joinpath(file_name)
    file_path.touch()
    with file_path.open('w') as output:
        output.write(text)
    return file_path.resolve()


def join_path(dir_path, name_to_join):
    return str(Path(dir_path).joinpath(name_to_join))


def change_file_extension(path, suffix):
    return str(Path(path).with_suffix(suffix))


def get_real_path(relativePath):
    return Path().cwd().joinpath(relativePath).resolve()


def copy_dir(source, destination):
    if Path(destination).exists():
        remove_dir(destination)
    # Path.mkdir(Path(destination))
    shutil.copytree(source, destination)

