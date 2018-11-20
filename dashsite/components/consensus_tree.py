import dash_core_components as dcc
import plotly.graph_objs as go

from igraph import *
import numpy as np
from pang.fileformat.json.JSONPangenome import JSONPangenome
import networkx as nx
from networkx.readwrite import json_graph


def get_tree(jsonpangenome: JSONPangenome, clck_data):
    G = nx.DiGraph()

    G.add_node(0)
    G.add_node(1)
    G.add_node(2)
    G.add_edge(0, 1)
    G.add_edge(0, 2)
    # G.add_node("ROOT")
    #
    # for i in range(5):
    #     G.add_node("Child_%i" % i)
    #     G.add_node("Grandchild_%i" % i)
    #     G.add_node("Greatgrandchild_%i" % i)
    #
    #     G.add_edge("ROOT", "Child_%i" % i)
    #     G.add_edge("Child_%i" % i, "Grandchild_%i" % i)
    #     G.add_edge("Grandchild_%i" % i, "Greatgrandchild_%i" % i)
    data = json_graph.tree_data(G, root=0)
    H = json_graph.tree_graph(data)
    return G



def get_consensus_tree_graph(jsonpangenome: JSONPangenome, networkx_tree):
    #todo clickdata type
    # create tree in networkx
    # if clickData is not None
    # read positions
    nr_vertices = 5
    v_label = map(str, range(nr_vertices))
    G = Graph.Tree(nr_vertices, 2)  # 2 stands for children number
    lay = G.layout('rt')

    position = {k: lay[k] for k in range(nr_vertices)}
    position = {0: [1,1], 1: [2,2], 2: [10,10], 3: [3, 3], 4: [3, 2]}
    Y = [lay[k][1] for k in range(nr_vertices)]
    M = max(Y)

    es = EdgeSeq(G)  # sequence of edges
    E = [e.tuple for e in G.es]  # list of edges

    L = len(position)
    Xn = [position[k][0] for k in range(L)]
    Yn = [2 * M - position[k][1] for k in range(L)]
    Xe = []
    Ye = []
    for edge in E:
        Xe += [position[edge[0]][0], position[edge[1]][0], None]
        Ye += [2 * M - position[edge[0]][1], 2 * M - position[edge[1]][1], None]

    labels = v_label

    lines = go.Scatter(x=Xe,
                       y=Ye,
                       mode='lines',
                       line=dict(color='rgb(210,210,210)', width=1),
                       hoverinfo='none'
                       )
    dots = go.Scatter(x=Xn,
                      y=Yn,
                      mode='markers',
                      name='',
                      marker=dict(symbol='circle',
                                  size=18,
                                  color='#6175c1',  # '#DB4551',
                                  line=dict(color='rgb(50,50,50)', width=1)
                                  ),
                      text=list(labels),
                      hoverinfo='text',
                      opacity=0.8
                      )
    # def make_annotations(pos, text, font_size=10, font_color='rgb(250,250,250)'):
    #     L = len(pos)
    #     if len(text) != L:
    #         raise ValueError('The lists pos and text must have the same len')
    #     annotations = go.layout.Annotation()
    #     for k in range(L):
    #         annotations.append(
    #             go.layout.Annotation(
    #                 text=labels[k],  # or replace labels with a different list for the text within the circle
    #                 x=pos[k][0], y=2 * M - position[k][1],
    #                 xref='x1', yref='y1',
    #                 font=dict(color=font_color, size=font_size),
    #                 showarrow=False)
    #         )
    #     return annotations

    axis = dict(showline=False,  # hide axis line, grid, ticklabels and  title
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                )

    layout = dict(title='Consensuses Tree',
                  # annotations=make_annotations(position, v_label),
                  font=dict(size=12),
                  showlegend=False,
                  xaxis=go.layout.XAxis(axis),
                  yaxis=go.layout.YAxis(axis),
                  margin=dict(l=40, r=40, b=85, t=100),
                  hovermode='closest',
                  plot_bgcolor='rgb(248,248,248)'
                  )

    # data = go.Data([lines, dots])
    # fig = dict(data=data, layout=layout)
    # fig['layout'].update(annotations=make_annotations(position, v_label))
    # return dcc.Graph(
    #     figure=go.Figure(
    #         data=data,
    #         layout=layout
    #     )
    # )
    # return go.Figure(
    #         data=[lines, dots],
    #         layout=layout
    #     )
    return {
            'data': [
                {
                    'x': [40, 20, 30, 40],
                    'y': [40, 1, 9, 50],
                    'text': ['aaa', 'b', 'c', 'd'],
                    'customdata': ['c.a', 'c.b', 'c.c', 'c.d'],
                    'name': 'Trace 1',
                    'mode': 'markers',
                    'marker': {'size': 12}
                },
                {
                    'x': [1, 2, 3, 4],
                    'y': [9, 4, 1, 4],
                    'text': ['w', 'x', 'y', 'z'],
                    'customdata': ['c.w', 'c.x', 'c.y', 'c.z'],
                    'name': 'Trace 2',
                    'mode': 'markers',
                    'marker': {'size': 12}
                }
            ]
        }


# def get_consensus_tree_scatter(jsonified_pangenome):
#     #simple scatter
#     N = 1000
#     random_x = np.random.randn(N)
#     random_y = np.random.randn(N)
#
#     # Create a trace
#     trace = go.Scatter(
#         x=random_x,
#         y=random_y,
#         mode='markers'
#     )
#     data = [trace]
#     return dcc.Graph(
#         figure=go.Figure(
#             data=data
#         )
#     )
#
#
# def get_consensus_tree1(jsonified_pangenome):
#
#     return dcc.Graph(
#         figure=go.Figure(
#             data=[
#                 go.Bar(
#                     x=[1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003,
#                        2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012],
#                     y=[219, 146, 112, 127, 124, 180, 236, 207, 236, 263,
#                        350, 430, 474, 526, 488, 537, 500, 439],
#                     name='Rest of world',
#                     marker=go.bar.Marker(
#                         color='rgb(55, 83, 109)'
#                     )
#                 ),
#                 go.Bar(
#                     x=[1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003,
#                        2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012],
#                     y=[16, 13, 10, 11, 28, 37, 43, 55, 56, 88, 105, 156, 270,
#                        299, 340, 403, 549, 499],
#                     name='China',
#                     marker=go.bar.Marker(
#                         color='rgb(26, 118, 255)'
#                     )
#                 )
#             ],
#             layout=go.Layout(
#                 title='US Export of Plastic Scrap',
#                 showlegend=True,
#                 legend=go.layout.Legend(
#                     x=0,
#                     y=1.0
#                 ),
#                 margin=go.layout.Margin(l=40, r=0, t=40, b=30)
#             )
#         ),
#         style={'height': 300},
#         id='my-graph'
#     )