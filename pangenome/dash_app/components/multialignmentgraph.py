from typing import List, Dict

from pangenome.pang.fileformats.json.JSONPangenome import JSONPangenome
import plotly.graph_objs as go


def get_data(jsonpangenome: JSONPangenome) -> Dict:
    multialignmentgraph_data = {}
    for sequence in jsonpangenome.sequences:
        multialignmentgraph_data[sequence.sequence_int_id] = {
            "seqid": sequence.sequence_str_id,
            "x": [jsonpangenome.nodes[node_id].column_id for node_id in sequence.nodes_ids],
            "y": [-sequence.sequence_int_id] * len(sequence.nodes_ids),
            "nodes_ids": [node_id for node_id in sequence.nodes_ids],
            "nodes_bases": [jsonpangenome.nodes[node_id].base for node_id in sequence.nodes_ids],
            "info": []
        }
    return multialignmentgraph_data


def get_graph(multialignmentgraph_data: Dict) -> go.Figure:

    sequences_traces = []
    sequences_names = []
    sequences_ids = []
    for seqid, sequencegraph_data in multialignmentgraph_data.items():
        sequences_ids.append(sequencegraph_data["y"][0])
        sequences_names.append(sequencegraph_data["seqid"])
        sequences_traces.extend(get_scattergls(x_pos=sequencegraph_data["x"],
                                               y_pos=sequencegraph_data["y"],
                                               customdata=sequencegraph_data["info"],
                                               labels=sequencegraph_data["nodes_bases"],
                                               hoverdata=sequencegraph_data["nodes_ids"]
                                               ))

    layout = dict(title='Multialignment Graph',
                  annotations=[],
                  font=dict(size=12),
                  showlegend=False,
                  xaxis=go.layout.XAxis(dict(title="Column ID", showline=False, zeroline=False, showgrid=False)),
                  yaxis=go.layout.YAxis(dict(title="SEQID", ticktext=sequences_names,
                                             tickvals=sequences_ids,
                                             showline=False, zeroline=False, showgrid=False, showticklabels=True,)),
                  margin=dict(l=40, r=40, b=85, t=100),
                  hovermode='closest',
                  plot_bgcolor='rgb(248,248,248)',
                  autosize=True,
                  )

    return go.Figure(
            data=sequences_traces,
            layout=layout
            )

def get_scattergls(x_pos: List[int], y_pos: List[int], customdata: List[str], labels: List[str], hoverdata: List[str]) -> List[go.Scattergl]:
    return [go.Scattergl(
        x=x_pos,
        y=y_pos,
        customdata=customdata,
        mode="markers",
        hoverinfo="text",
        text=hoverdata,
        # marker="circle",
        marker=dict(
            symbol="circle",
            size=15,
            color='white',
            line=dict(
                width=1,
                color='rgb(0,0,0)'),
        ),
        textfont=dict(
            size=12,
            color='rgb(255,255,255)'
        )
    )
    , go.Scattergl(
        x=x_pos,
        y=y_pos,
        customdata=customdata,
        mode='text',
        text=labels,
        textfont=dict(
            size=12
        )
    )]

