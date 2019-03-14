from dash.dependencies import Input, Output, State

from dash_app.components import consensustable, consensustree
from ..components import jsontools
from ..layout.layout_ids import *

from ..server import app

@app.callback(
    Output(id_partial_consensustable_hidden, 'children'),
    [Input(id_full_consensustable_hidden, 'children'),
     Input(id_full_consensustree_hidden, 'children'),
     Input(id_consensus_tree_slider, 'value')]
)
def update_partial_table_data(jsonified_full_consensustable: str, jsonified_tree: str, slider_value: float):
    if not jsonified_full_consensustable or not jsonified_tree:
        return []
    full_consensustable_data = jsontools.unjsonify_df(jsonified_full_consensustable)
    full_consensustree_data = jsontools.unjsonify_dict(jsonified_tree)
    full_consensustree_tree = consensustree.dict_to_tree(full_consensustree_data)
    table_without_consensuses_smaller_than_slider = consensustable.remove_smaller_than_slider(full_consensustable_data,
                                                                                              full_consensustree_tree,
                                                                                              slider_value)
    # tu dołożyć ewentualne zmiany wynikające z innych elementów Input/State
    return jsontools.jsonify_df(table_without_consensuses_smaller_than_slider)

@app.callback(
    Output(id_consensuses_table, 'data'),
    [Input(id_partial_consensustable_hidden, 'children')]
)
def to_consensustable_content(jsonified_partial_consensustable):
    if not jsonified_partial_consensustable:
        return []
    partial_consensustable_data = jsontools.unjsonify_df(jsonified_partial_consensustable)
    data_rows = partial_consensustable_data.to_dict("rows")
    return data_rows

@app.callback(
    Output(id_consensuses_table, 'columns'),
    [Input(id_partial_consensustable_hidden, 'children')]
)
def update_columns(jsonified_partial_consensustable):
    if not jsonified_partial_consensustable:
        return [{}]
    partial_consensustable_data = jsontools.unjsonify_df(jsonified_partial_consensustable)
    return [{"name": i, "id": i} for i in partial_consensustable_data.columns]

@app.callback(
    Output(id_consensuses_table, 'style_data_conditional'),
    [Input(id_partial_consensustable_hidden, 'children')],
    [State(id_full_consensustree_hidden, 'children')]
)
def color_consensuses_table_cells(jsonified_partial_consensustable, jsonified_consensus_tree):
    if not jsonified_partial_consensustable or not jsonified_consensus_tree:
        return []
    partial_consensustable_data = jsontools.unjsonify_df(jsonified_partial_consensustable)
    consensustree_data = jsontools.unjsonify_dict(jsonified_consensus_tree)
    tree = consensustree.dict_to_tree(consensustree_data)

    return consensustable.get_cells_styling(tree, partial_consensustable_data)
