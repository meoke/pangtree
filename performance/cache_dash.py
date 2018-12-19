import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input
from dash.dependencies import Output
import numpy as np
# import json
import jsonpickle

# global_df = pd.read_csv('...')
app = dash.Dash(__name__)
app.layout = html.Div([
    html.Div(id='a'),
    dcc.Input(id='i'),
    html.Div(id='h')
])

@app.callback(Output('a', 'children'),
              [Input('h', 'children')])
def show_in_div(jsonified_value):
    a_v  = jsonpickle.decode(jsonified_value)
    a = [a_v, a_v]
    return f"from hidden: {a}"


@app.callback(Output('h', 'children'),
              [Input('i', 'value')])
def show_in_div(value):
     return jsonpickle.encode(int(value))


if __name__ == '__main__':
    app.run_server(debug=True, port=1080)
