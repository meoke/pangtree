import pandas as pd
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input
from dash.dependencies import Output
import dash_table
import numpy as np
# import json
import jsonpickle

# global_df = pd.read_csv('...')
app = dash.Dash(__name__)
d = {'Region': [1, 2], 'Continent': [3, 4]}
df = pd.DataFrame(data=d)
app.layout = html.Div([
    html.Div(id='a'),
    dcc.Input(id='i'),
    html.Div(id='h'),
    dash_table.DataTable(
        data=df.to_dict('rows'),
        columns=[
            {'name': i, 'id': i} for i in df.columns
        ],
        style_data_conditional=[
            {
                'if': {
                    'column_id': 'Region',
                    'filter': 'Region eq "Montreal"'
                },
                'backgroundColor': '#3D9970',
                'color': 'white',
            },
            {
                'if': {
                    'column_id': 'Humidity',
                    'filter': 'Humidity eq num(20)'
                },
                'backgroundColor': '#3D9970',
                'color': 'white',
            },
            {
                'if': {
                    'column_id': 'Temperature',
                    'filter': 'Temperature > num(3.9)'
                },
                'backgroundColor': '#3D9970',
                'color': 'white',
            },
        ]
    )
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
