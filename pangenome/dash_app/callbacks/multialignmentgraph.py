from dash.dependencies import Input, Output

from dash_app.components import multialignmentgraph
from ..components import jsontools
from ..layout.layout_ids import *

from ..server import app

@app.callback(
    Output(id_multialignmentgraph, 'figure'),
    [Input(id_multialignmentgraph_hidden, 'children')]
)
def update_multialignmentgraph(jsonified_multialignmentgraph_data: str):
    if not jsonified_multialignmentgraph_data:
        return []
    multialignmentgraph_data = jsontools.unjsonify_dict(jsonified_multialignmentgraph_data)
    return multialignmentgraph.get_graph(multialignmentgraph_data)


@app.callback(
    Output(id_multialignmentgraph_container, 'style'),
    [Input(id_multialignmentgraph_hidden, 'children')])
def show_multialignemntgraph(jsonified_multialignmentgraph_data):
    if jsonified_multialignmentgraph_data:
        return {'display': 'block'}
    else:
        return {'display': 'none'}
