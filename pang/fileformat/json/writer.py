from pathlib import Path
import Pangenome
from fileformat.json.JSONPangenome import JSONPangenome
from userio import pathtools
import jsonpickle


def save(output_dir: Path, pangenome: Pangenome):
    jsonpoagraph = JSONPangenome(pangenome)
    json_path = pathtools.get_child_file_path(output_dir, "pangenome.json")
    jsonpickle.set_encoder_options('simplejson', indent=4)
    with open(json_path, 'w') as json_output:
        json_output.write(jsonpickle.encode(jsonpoagraph, unpicklable=False))
    return json_path
