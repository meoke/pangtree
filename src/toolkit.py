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
            last_dir_suffix = int(existing_prefixed_dirs[-1].stem.split("_")[-1])
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

def get_next_child_file_name(path, child_file_prefix):
    existing_prefixed_files = sorted(list(path.glob("*"+child_file_prefix.split('.')[0]+"*")))
    if existing_prefixed_files:
        try:
            last_path_suffix = int(existing_prefixed_files[-1].stem.split("_")[-1])
        except:
            last_path_suffix = -1
    else:
        last_path_suffix = -1
    new_file_suffix_value = last_path_suffix + 1
    new_file_suffix = str(new_file_suffix_value) if new_file_suffix_value > 9 else '0' + str(new_file_suffix_value)
    new_file_name = child_file_prefix.split('.')[0] + '_' + new_file_suffix + "." + child_file_prefix.split('.')[1]
    new_file_path = path.joinpath(new_file_name)
    return str(new_file_path)

def create_child_dir(parent_path, dir_name):
    child_dir_path = Path(parent_path).joinpath(dir_name)
    if child_dir_path.exists():
        return create_next_sibling_dir(child_dir_path, dir_name)
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


def get_real_path(relative_path):
    return Path().cwd().joinpath(relative_path).resolve()


def copy_dir(source, destination):
    if Path(destination).exists():
        remove_dir(destination)
    shutil.copytree(source, destination)


def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)