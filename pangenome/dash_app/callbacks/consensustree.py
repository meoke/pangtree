from dash.dependencies import Input, Output, State
import dash_html_components as html

from dash_app.components import consensustable
from ..components import jsontools
from ..layout.layout_ids import *
from ..components import consensustree
from ..server import app

@app.callback(
    Output(id_current_consensustree_hidden, 'children'),
    [Input(id_full_consensustree_hidden, 'children')]
)
def update_current_tree_state(jsonified_full_consensustree):
    if not jsonified_full_consensustree:
        return []
    full_consensustree_data = jsontools.unjsonify_dict(jsonified_full_consensustree)
    full_consensus_tree = consensustree.dict_to_tree(full_consensustree_data)
    # any changes to the tree put here

    current_consensustree_data = consensustree.tree_to_dict(full_consensus_tree)
    return jsontools.jsonify_dict(current_consensustree_data)


@app.callback(
    Output(id_consensus_tree_graph, 'figure'),
    [Input(id_current_consensustree_hidden, 'children'),
     Input(id_consensus_tree_slider, 'value'),
     Input(id_leaf_info_dropdown, 'value'),
     Input(id_full_consensustable_hidden, 'children')])
def to_consensustree_graph(jsonified_current_consensustree, slider_value, leaf_info, jsonified_full_consensustable):
    if not jsonified_current_consensustree or not jsonified_full_consensustable:
        return []
    current_consensustree_data = jsontools.unjsonify_dict(jsonified_current_consensustree)
    current_consensustree_tree = consensustree.dict_to_tree(current_consensustree_data)
    full_consensustable_data = jsontools.unjsonify_df(jsonified_full_consensustable)
    graph = consensustree.get_consensustree_graph(current_consensustree_tree, slider_value, leaf_info, full_consensustable_data)
    return graph

@app.callback(
    Output(id_consensus_tree_container, 'style'),
    [Input(id_current_consensustree_hidden, 'children')])
def show_consensus_tree_container(jsonified_current_consensustree):
    if jsonified_current_consensustree:
        return {'display': 'block'}
    else:
        return {'display': 'none'}

@app.callback(
    Output(id_leaf_info_dropdown, 'options'),
    [Input(id_full_consensustable_hidden, 'children')])
def to_consensustree_leaf_info_options_dropdown(jsonified_full_consensustable):
    if not jsonified_full_consensustable:
        return []
    full_consensustable = jsontools.unjsonify_df(jsonified_full_consensustable)
    metadata = consensustable.get_metadata_list(full_consensustable)
    dropdown_options = consensustree.get_leaf_info_dropdown_options(metadata)
    return dropdown_options