from dash.dependencies import Input, Output

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
    Output(id_consensus_node_details_header, 'children'),
    [Input(id_consensus_tree_graph, 'clickData')]
)
def to_consensus_node_details_header(tree_click_data):
    if not tree_click_data:
        return []
    clicked_node = tree_click_data['points'][0]
    node_id = clicked_node['pointIndex']
    return f"Consensus {node_id}"

@app.callback(
    Output(id_consensus_node_details_table_hidden, 'children'),
    [Input(id_consensus_tree_graph, 'clickData'),
     Input(id_full_consensustable_hidden, 'children'),
     Input(id_full_consensustree_hidden, 'children')]
)
def to_consensus_node_details_table(tree_click_data, jsonified_full_consensustable, jsonified_consensustree):
    if not jsonified_full_consensustable or not tree_click_data:
        return []
    clicked_node = tree_click_data['points'][0]
    node_id = clicked_node['pointIndex']
    full_consensustable = jsontools.unjsonify_df(jsonified_full_consensustable)
    consensustree_data = jsontools.unjsonify_dict(jsonified_consensustree)
    tree = consensustree.dict_to_tree(consensustree_data)
    node_details_df = consensustable.get_consensus_details_df(node_id, full_consensustable, tree)
    return jsontools.jsonify_df(node_details_df)

@app.callback(
    Output(id_consensus_node_details_table, 'data'),
    [Input(id_consensus_node_details_table_hidden, 'children')]
)
def to_consensusnode_details_content(jsonified_consensus_details_table):
    if not jsonified_consensus_details_table:
        return []
    consensus_details_table_data = jsontools.unjsonify_df(jsonified_consensus_details_table)
    data_rows = consensus_details_table_data.to_dict("rows")
    return data_rows

@app.callback(
    Output(id_consensus_node_details_table, 'columns'),
    [Input(id_consensus_node_details_table_hidden, 'children')]
)
def update_columns(jsonified_consensus_details_table):
    if not jsonified_consensus_details_table:
        return [{}]
    consensus_details_table_data = jsontools.unjsonify_df(jsonified_consensus_details_table)
    return [{"name": i, "id": i} for i in list(consensus_details_table_data.columns)]

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