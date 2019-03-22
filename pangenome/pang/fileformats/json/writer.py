from pathlib import Path
from Pangenome import Pangenome
from .JSONPangenome import JSONPangenome
from tools import pathtools
import jsonpickle

#
# def pangenome_to_jsonpangenome(pangenome: Pangenome):
#     return JSONPangenome(pangenome)
#
#
# def pangenome_to_json(pangenome: Pangenome):
#     jsonpangenome = pangenome_to_jsonpangenome(pangenome)
#     return jsonpickle.encode(jsonpangenome)


def save_to_file(output_dir: Path, pangenome: Pangenome):
    jsonpoagraph = JSONPangenome(pangenome, pangenome.params)
    json_path = pathtools.get_child_path(output_dir, "pangenome.json")
    jsonpickle.set_encoder_options('simplejson', indent=4)
    with open(json_path, 'w') as json_output:
        json_output.write(jsonpickle.encode(jsonpoagraph, unpicklable=True))
    return json_path

#
# def jsonpangenome_to_json(jsonpangenome: JSONPangenome):
#     return jsonpickle.encode(jsonpangenome)
