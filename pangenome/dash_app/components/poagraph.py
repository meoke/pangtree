from typing import List, Dict, Union

from fileformats.json.JSONPangenome import JSONPangenome
import plotly.graph_objs as go

y_pos_dict = {'A': 40, 'C': 35, 'G': 30, 'T': 25, 'N': 20, 'W': 10, '?': 5, 'n': 2}

PointsDict = Dict[str, Union[List[int], List[str]]]
Path = Dict[str, Union[List[int], List[str]]]
PathsDict = Dict[int, Path]
PoagraphData = Dict[str, Union[PointsDict, PathsDict]]


def get_data(jsonpangenome: JSONPangenome) -> PoagraphData:
    poagraph_data: PoagraphData = {"nodes": {"x": [node.column_id for node in jsonpangenome.nodes],
                                             "y": [y_pos_dict[node.nucleobase] for node in jsonpangenome.nodes],
                                             "base": [node.nucleobase for node in jsonpangenome.nodes]},
                                   "paths": {}}

    for sequence in jsonpangenome.sequences:
        poagraph_data["paths"][sequence.sequence_int_id] = {
            "x": [jsonpangenome.nodes[node_id].column_id for node_id in sequence.nodes_ids],
            "y": [y_pos_dict[jsonpangenome.nodes[node_id].nucleobase] for node_id in sequence.nodes_ids],
            "name": sequence.sequence_str_id
        }

    return poagraph_data


def get_graph(poagraph_data: Dict) -> go.Figure:
    nodes_scatter = get_nodes_scatter(poagraph_data["nodes"])
    scatters = []

    for seqid, sequencegraph_data in poagraph_data["paths"].items():
        scatters.append(get_path_scatter(name=sequencegraph_data["name"],
                                         x_pos=sequencegraph_data["x"],
                                         y_pos=sequencegraph_data["y"]
                                         ))
    scatters.append(nodes_scatter)
    layout = dict(title='Poagraph',
                  annotations=[],
                  font=dict(size=12),
                  showlegend=True,
                  xaxis=go.layout.XAxis(dict(title="Column ID", showline=False, zeroline=False, showgrid=False)),
                  yaxis=go.layout.YAxis(dict(title="Base",
                                             ticktext=[*y_pos_dict.keys()],
                                             tickvals=[*y_pos_dict.values()],
                                             showline=False, zeroline=False, showgrid=False, showticklabels=True,)),
                  margin=dict(l=40, r=40, b=85, t=100),
                  hovermode='closest',
                  plot_bgcolor='rgb(248,248,248)',
                  autosize=True,
                  )

    return go.Figure(
            data=scatters,
            layout=layout
            )

def get_nodes_scatter(nodes_poagraph_data) -> go.Scatter:
    return go.Scatter(x=nodes_poagraph_data["x"],
                      y=nodes_poagraph_data["y"],
                      text=nodes_poagraph_data["base"],
                      mode="markers+text",
                      marker=dict(
                          symbol="circle",
                          size=20,
                          color='white',
                          line=dict(
                              width=1),
                      ),
                      textfont=dict(
                          size=12,
                          color='rgb(0,0,0)'
                      )
                      )


def get_path_scatter(name: str, x_pos: List[str], y_pos: List[str]) -> go.Scatter:
    return go.Scatter(
        x=x_pos,
        y=y_pos,
        name=name,
        mode='lines',
            line=dict(
                    width=1,
                    color='#404040',
                    shape='spline'
            )
    )


