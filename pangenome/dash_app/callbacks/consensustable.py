from dash.dependencies import Input, Output
import dash_html_components as html
from ..components import jsontools
from ..layout.layout_ids import *

from ..server import app

@app.callback(
    Output(id_partial_consensustable_hidden, 'children'),
    [Input(id_full_consensustable_hidden, 'children')]
)
def update_partial_table_data(jsonified_full_consensustable):
    if not jsonified_full_consensustable:
        return []
    full_consensustable_data = jsontools.unjsonify_df(jsonified_full_consensustable)
    # tu dołożyć ewentualne zmiany wynikające z innych elementów Input/State
    return jsontools.jsonify_df(full_consensustable_data)

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
        return []
    partial_consensustable_data = jsontools.unjsonify_df(jsonified_partial_consensustable)
    return [{"name": i, "id": i} for i in partial_consensustable_data.columns]