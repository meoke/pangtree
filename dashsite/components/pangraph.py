import pandas as pd
import plotly.graph_objs as go


def get_graph(nodes_data, traces_data):
    g = {
            'data': [go.Scattergl(
                x=nodes_data.x,
                y=nodes_data.y,
                mode='markers',
                hoverinfo='text',
                text=nodes_data.ID,
                # marker=dict(
                #     line=dict(
                #         width=1,
                #         color='#404040')
                # ),
                # line=dict(
                #     width=[1,5,10]
                # ),
                name='nodes',
                # ids=["1","2",30]
                )] +
    [ go.Scattergl(
        x=seq_nodes[0],
        y=seq_nodes[1],
        mode='lines',
        # line = dict(
        #     shape='spline'
        # ),
        hovertext="",
        hoverinfo='skip',
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
                'title': 'Pangraph',
                # 'width': nodes_data.shape[0]/2,
                'hovermode': 'closest'
            }
        }
    print("scatters ready")
    return g

def get_traces_data(jsonpangenome):
    y_position_dict = {'A': 40, 'C': 35, 'G': 30, 'T': 25, 'N': 20, 'W': 10}
    d = {
        sequence.mafname: ([jsonpangenome.nodes[node_id].column_id * 5 for node_id in sequence.nodes_ids],
                            [y_position_dict[jsonpangenome.nodes[node_id].nucleobase] for node_id in sequence.nodes_ids])
        for sequence in jsonpangenome.sequences
    }
    return d

def get_nodes_data(jsonpangenome):
    y_position_dict = {'A': 40, 'C': 35, 'G': 30, 'T': 25, 'N': 20, 'W': 10}

    nodes_data = pd.DataFrame(
        columns=['ID', 'x', 'y', 'base', 'block_id', 'aligned_to'])

    nodes_data_list = [
        {"ID": node.id,
         "x": node.column_id * 5,
         "y": y_position_dict[node.nucleobase],
         "base": node.nucleobase,
         "block_id": node.block_id,
         "aligned_to": node.aligned_to}
        for node in jsonpangenome.nodes
    ]
    nodes_data = nodes_data.append(nodes_data_list, ignore_index=True)
    print("nodes data ready")
    return nodes_data
