from dash.dependencies import Input, Output, State

from dash_app.components import parameters, consensustable

import dash_app.components.jsontools as jsontools
from ..layout.layout_ids import *

from ..server import app

@app.callback(
    Output(id_pangenome_parameters_hidden, 'children'),
    [Input(id_pangenome_hidden, 'children')]
)
def update_pangenome_parameters_hidden(jsonified_pangenome):
    if not jsonified_pangenome:
        return ""
    jsonpangenome = jsontools.unjsonify_jsonpangenome(jsonified_pangenome)
    parameters_data = parameters.get_data(jsonpangenome)
    return jsontools.jsonify_dict(parameters_data)


@app.callback(
    Output(id_full_consensustable_hidden, 'children'),
    [Input(id_pangenome_hidden, 'children')]
)
def update_full_consensustable_hidden(jsonified_pangenome):
    if not jsonified_pangenome:
        return []
    jsonpangenome = jsontools.unjsonify_jsonpangenome(jsonified_pangenome)
    consensustable_data = consensustable.get_full_table_data(jsonpangenome)
    return jsontools.jsonify_df(consensustable_data)


@app.callback(
    Output(id_pangenome_consensustree_hidden, 'children'),
    [Input(id_pangenome_hidden, 'children')]
)
def update_consensustree_hidden(jsonified_pangenome):
    return "Pełne dane o consensus tree"

@app.callback(
    Output(id_pangraph_hidden, 'children'),
    [Input(id_pangenome_hidden, 'children')]
)
def update_pangraph_hidden(jsonified_pangenome):
    return "Pełne dane o pangraph"

@app.callback(
    Output(id_mafgraph_hidden, 'children'),
    [Input(id_pangenome_hidden, 'children')]
)
def update_mafgraph_hidden(jsonified_pangenome):
    return "Pełne dane o mafgraph"