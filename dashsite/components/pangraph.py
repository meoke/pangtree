import pandas as pd
import plotly.graph_objs as go
import numpy as np


def get_graph(pangraph_data):
    return {
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
                    x=np.array(range(10)),
                    y=np.sin(np.array(range(10))),
                    mode='lines',
                    line=dict(
                            width=1,
                            color='#404040',
                            shape='spline'
                            ),
                    opacity=0.2,
                    name='DWA'
                )
            ],
            'layout': {
                'title': 'Pangraph'
            }
        }

def get_data(jsonified_pangenome):
    d = {'col1': [1, 2], 'col2': [3, 4]}
    df = pd.DataFrame(data=d)
    return df