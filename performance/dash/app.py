import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

N = 1000000

import dash
import dash_core_components as dcc
import dash_html_components as html

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dcc.Graph(
        id='example-graph',
        figure={
            'data': [go.Scattergl(
                x=[0,-1,2],
                y=[0,1,2],
                mode='lines',
                marker=dict(
                    line=dict(
                        width=1,
                        color='#404040')
                ),
                name='JEDEN',
                ids=["1","2",30]
            ),
                go.Scatter(
                    x=np.array(range(100000)),
                    y=np.sin(np.array(range(100000))),
                    mode='lines',
                    line=dict(
                            width=1,
                            color='#404040',
                            shape='spline'
                            ),
                    opacity=0.2,
                    # marker=dict(
                    #     line=dict(
                    #         width=10,
                    #         color='#404040',
                    #         shape='hv'
                    #         )

                    name='DWA'
                ),
                go.Scattergl(
                    x=[0, 1, 2],
                    y=[0, 1, 2],
                    mode='lines',
                    marker=dict(
                        line=dict(
                            width=1,
                            color='#404040')
                    ),
                    name='TRZY'
                )
            ],
            'layout': {
                'title': 'Dash Data Visualization'
            }
        }
    )
])
app.scripts.config.serve_locally=True

if __name__ == '__main__':
    app.run_server(debug=True)
