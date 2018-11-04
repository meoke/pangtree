from pathlib import Path
from graph import Pangraph
from metadata import MultialignmentMetadata
from userio import pathtools
from .JSONPoagraph import JSONPangraph
import jsonpickle


def save(output_dir: Path, pangraph: Pangraph, genomes_info: MultialignmentMetadata):
    jsonpoagraph = JSONPangraph(pangraph, genomes_info)
    json_path = pathtools.get_child_file_path(output_dir, "pangraph.json")
    with open(json_path, 'w') as json_output:
        json_output.write(jsonpickle.dumps(jsonpoagraph, unpicklable=False))
