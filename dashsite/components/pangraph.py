import pandas as pd
import plotly.graph_objs as go
import numpy as np


def get_graph(nodes_data, traces_data):
    g = {
            'data': [go.Scattergl(
                x=nodes_data.x,
                y=nodes_data.y,
                mode='markers',
                # marker=dict(
                #     line=dict(
                #         width=1,
                #         color='#404040')
                # ),
                name='nodes',
                # ids=["1","2",30]
                )] +
    [ go.Scattergl(
        x=list(zip(*seq_nodes))[0],
        y=list(zip(*seq_nodes))[1],
        mode='lines',
        # line = dict(
        #     shape='spline'
        # ),
        name=seq_name
    )
                # go.Scatter(
                #     x=np.array(range(10)),
                #     y=np.sin(np.array(range(10))    ),
                #     mode='lines',
                #     line=dict(
                #             width=1,
                #             color='#404040',
                #             shape='spline'
                #             ),
                #     opacity=0.2,
                #     name='DWA'
                # )
             for seq_name, seq_nodes in traces_data.items()],
            'layout': {
                'title': 'Pangraph'
            }
        }
    print("scatters ready")
    return g

def get_traces_data(jsonpangenome):
    y_position_dict = {'A': 40, 'C': 30, 'G': 20, 'T': 10, 'N': 0}
    d = {
        sequence.mafname: [(jsonpangenome.nodes[node_id].column_id * 10,
                            y_position_dict[jsonpangenome.nodes[node_id].nucleobase])
                           for node_id in sequence.nodes_ids]
        for sequence in jsonpangenome.sequences
    }
    return d

def get_nodes_data(jsonpangenome):
    y_position_dict = {'A': 40, 'C': 30, 'G': 20, 'T': 10, 'N': 0}

    nodes_data = pd.DataFrame(
        columns=['ID', 'x', 'y', 'base', 'block_id', 'aligned_to'])

    for node in jsonpangenome.nodes:
        node = {"ID": node.id,
                "x": node.column_id * 10,
                "y": y_position_dict[node.nucleobase],
                "base": node.nucleobase,
                "block_id": node.block_id,
                "aligned_to": node.aligned_to}

        nodes_data = nodes_data.append(node, ignore_index=True)
    return nodes_data
