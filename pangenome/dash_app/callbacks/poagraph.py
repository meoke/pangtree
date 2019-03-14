from dash.dependencies import Input, Output

from dash_app.components import multialignmentgraph, poagraph
from ..components import jsontools
from ..layout.layout_ids import *

from ..server import app

@app.callback(
    Output(id_poagraph, 'figure'),
    [Input(id_poagraph_hidden, 'children')]
)
def update_poagraph(jsonified_poagraph_data: str):
    if not jsonified_poagraph_data:
        return []
    poagraph_data = jsontools.unjsonify_dict(jsonified_poagraph_data)
    return poagraph.get_graph(poagraph_data)


@app.callback(
    Output(id_poagraph_container, 'style'),
    [Input(id_poagraph_hidden, 'children')])
def show_poagraph(jsonified_poagraph_data):
    if jsonified_poagraph_data:
        return {'display': 'block'}
    else:
        return {'display': 'none'}
