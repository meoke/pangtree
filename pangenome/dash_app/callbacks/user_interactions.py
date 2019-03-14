import shutil
from base64 import b64decode
from pathlib import Path

import jsonpickle

from Pangenome import Pangenome
from arguments.PangenomeParameters import PangenomeParameters
from fileformats.json.JSONPangenome import JSONPangenome
from ..server import app
from dash.dependencies import Input, Output, State
from ..layout.layout_ids import *
import json

@app.callback(
    Output(id_last_clicked_hidden, 'children'),
    [Input(id_pang_button, 'n_clicks'),
     Input(id_load_pangenome_button, 'n_clicks')],
    [State(id_last_clicked_hidden, 'children')]
)
def update_last_clicked_info(pang_n_clicks: int, load_pangenome_n_clicks: int, last_clicked_jsonified: str) -> str:
    new_clicked_info = {"pang": 0, "load": 0, "action": ""}
    if last_clicked_jsonified is None:
        if pang_n_clicks == 1:
            new_clicked_info = {"pang": 1, "load": 0, "action": "pang"}
        if load_pangenome_n_clicks == 1:
            new_clicked_info = {"pang": 0, "load": 1, "action": "load"}
    else:
        last_clicked_info = json.loads(last_clicked_jsonified)
        if pang_n_clicks and last_clicked_info["pang"] == pang_n_clicks - 1:
            new_clicked_info = {"pang": last_clicked_info["pang"]+1,
                               "load": last_clicked_info["load"],
                               "action": "pang"}
        if load_pangenome_n_clicks and last_clicked_info["load"] == load_pangenome_n_clicks - 1:
            new_clicked_info = {"pang": last_clicked_info["pang"],
                               "load": last_clicked_info["load"]+1,
                               "action": "load"}
    return json.dumps(new_clicked_info)

@app.callback(
    Output(id_pangenome_hidden, 'children'),
    [Input(id_last_clicked_hidden, 'children')],
    [State('pangenome_upload', 'contents'),
     State('maf_upload', 'contents'),
     State('metadata_upload', 'contents'),
     State('pang_options', 'values'),
     State('algorithm', 'value'),
     State('hbmin', 'value'),
     State('r', 'value'),
     State('multiplier', 'value'),
     State('stop', 'value'),
     State('tree_consensus_options', 'values')
     ])
def call_pang(last_clicked_jsonified: str,
              pangenome_contents,
              maf_contents,
              metadata_contents,
              pang_options_values,
              consensus_algorithm_value,
              hbmin_value,
              r_value,
              multiplier_value,
              stop_value,
              tree_consensus_options_values) -> str:
    last_clicked = json.loads(last_clicked_jsonified)
    if last_clicked["action"] == 'load':
       return decode_content(pangenome_contents)
        # return get_jsonified_pangenome_from_jsonpangenomefile(pangenome_decoded)
    elif last_clicked["action"] == "pang":
        fasta_option = True if 'FASTA' in pang_options_values else False
        re_consensus_value = True if 're_consensus' in tree_consensus_options_values else False
        anti_fragmentation_value = True if 'anti_fragmentation' in tree_consensus_options_values else False
        no_multiplier_anti_granular = True  # todo
        output_path = "" #todo
        maf_decoded = decode_content(maf_contents)
        metadata_decoded = decode_content(metadata_contents) if metadata_contents else None
        params = PangenomeParameters(multialignment_file_content=maf_decoded,
                                     multialignment_file_path="",
                                     metadata_file_content=metadata_decoded,
                                     metadata_file_path="",
                                     blosum_file_path="",
                                     output_path=output_path,
                                     generate_fasta=fasta_option,
                                     consensus_type="",
                                     hbmin=hbmin_value,
                                     search_range=r_value,
                                     multiplier=multiplier_value,
                                     stop=stop_value,
                                     re_consensus=re_consensus_value,
                                     not_dag="",
                                     fasta_complementation_option="",
                                     missing_nucleotide_symbol="",
                                     local_fasta_dirpath="",
                                     max_cutoff_option="",
                                     node_cutoff_option="")
        pangenome = run_pangenome_algorithm(params)
        jsonpangenome = get_jsonpangenome_from_pangenome(pangenome, params)
        jsonified_pangenome = get_jsonified_pangenome_from_jsonpangenome(jsonpangenome)
        pangenome_json_path = save_jsonifiedpangenome_to_file(jsonified_pangenome)
        shutil.copy(pangenome_json_path, "download/pangenome.json")
        return jsonified_pangenome
    else:
        return ""

def get_jsonified_pangenome_from_jsonpangenomefile(jsonpangenomefile_stream) -> str:
    content_type, content_string = jsonpangenomefile_stream.split(',')
    jsonified_pangenome = b64decode(content_string).decode('ascii')
    return jsonified_pangenome

def run_pangenome_algorithm(pangenome_parameters: PangenomeParameters) -> str:
    raise NotImplementedError("Run pangenome algorithm")

def get_jsonpangenome_from_pangenome(pangenome: Pangenome, pangenome_parameters: PangenomeParameters) -> JSONPangenome:
    return JSONPangenome(pangenome, pangenome_parameters)

def save_jsonifiedpangenome_to_file(jsonified_pangenome: str) -> Path:
    raise NotImplementedError("save_jsonifiedpangenome_to_file")
    jsonpath = Path("") #todo
    with open(jsonpath, 'w') as json_output:
        json_output.write(jsonified_pangenome)
    return jsonpath

def get_jsonified_pangenome_from_jsonpangenome(jsonpangenome: JSONPangenome) -> str:
    jsonpickle.set_encoder_option('simplejson', indent=4)
    return jsonpickle.encode(jsonpangenome, unpicklable=True)

def decode_content(content: str) -> str:
    if not content:
        return ''
    content_string = content.split(',')[1]
    return b64decode(content_string).decode('ascii')