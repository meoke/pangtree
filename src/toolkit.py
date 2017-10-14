#import os
from pathlib import Path
import shutil
# import itertools, collections
#
# def create_child_dir(currentDir, childDirName):
#     dirPath = os.path.join(currentDir, childDirName)
#     try:
#         os.mkdir(dirPath)
#         return(dirPath)
#     except:
#         return(dirPath)
#
#
# def getRealPath(relativePath):
#     return os.path.realpath(relativePath)
#
#
# def mergeFiles(inputDirPath, outputPath):
#     with open(outputPath, "w") as outfile:
#         for fname in os.listdir(inputDirPath):
#             if os.path.isfile(os.path.join(inputDirPath, fname)):
#                 filePath = os.path.join(inputDirPath, fname)
#                 with open(filePath) as infile:
#                     outfile.write(infile.read())
#
#


#
# def changeExtension(filePath, newExtension):
#     return os.path.splitext(filePath)[0] + newExtension
#
#
# def getBinPath(currentPath, programName):
#     return(os.path.join(currentPath, os.pardir, 'bin', programName))
#
#
# def consume(iterator, n):
#     "Advance the iterator n-steps ahead. If n is none, consume entirely."
#     # Use functions that consume iterators at C speed.
#     if n is None:
#         # feed the entire iterator into a zero-length deque
#         collections.deque(iterator, maxlen=0)
#     else:
#         # advance to the empty slice starting at position n
#         next(itertools.islice(iterator, n, n), None)
#
#

#
#
# def get_file_path(dir_path, file_name):
#     return os.path.join(dir_path, file_name)
#
#
# def save_file(file_path, file_content):
#     with open(file_path, 'w') as output:
#         output.write(file_content)
#

# def get_dir_name(dir_path):
#     return dir_path.split('/')[-1]
#
#
# def get_next_child_dir_name(full_ebola_file_path, dir_prefix):
#     if isdir(full_ebola_file_path):
#         current_dir_path = full_ebola_file_path
#     else:
#         current_dir_path = get_file_dirname(full_ebola_file_path)
#     child_dir_list = os.listdir(current_dir_path)
#     prefixed_dirs = sorted([*filter(lambda dir_name: True if dir_prefix in dir_name else  False, child_dir_list)])
#     if prefixed_dirs:
#         last_dir_number = int(prefixed_dirs[-1].split("_")[-1])
#     else:
#         last_dir_number = -1
#     if last_dir_number < 10:
#         next_dir_number = str(0) + str(last_dir_number+1)
#     else:
#         next_dir_number = str(last_dir_number + 1)
#     new_dir_name = "_".join([dir_prefix, next_dir_number])
#     create_child_dir(current_dir_path, new_dir_name)
#     return new_dir_name

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