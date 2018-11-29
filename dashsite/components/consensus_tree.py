import plotly.graph_objs as go

from pang.fileformat.json.JSONPangenome import JSONPangenome
import networkx as nx
from collections import deque
from networkx.drawing.nx_agraph import graphviz_layout


def toggle_node_with_children(old_tree, clicked_node):
    nodes_to_toggle = deque([clicked_node])
    while nodes_to_toggle:
        node_id = nodes_to_toggle.pop()
        old_tree.nodes[node_id]['hidden'] = not old_tree.nodes[node_id]['hidden']
        node_children = [v for u, v in old_tree.edges if u == node_id]
        nodes_to_toggle.extend(node_children)

    return old_tree


def get_tree(jsonpangenome: JSONPangenome, clck_data, old_tree):
    if not old_tree:
        tree = create_tree(jsonpangenome)
    else:
        tree = old_tree
    if clck_data:
        clicked_node = int(clck_data['points'][0]['customdata'])
        return toggle_node_with_children(tree, clicked_node)
    return tree


def create_tree(jsonpangenome: JSONPangenome):
    tree_graph = nx.DiGraph()
    for consensus in jsonpangenome.consensus_tree:
        tree_graph.add_node(consensus.id,
                            name=consensus.name,
                            comp=consensus.comp_to_all_sequences,
                            sequences=consensus.sequences,
                            show_in_table=True,
                            hidden=False,
                            children_consensuses=consensus.children,
                            mincomp=min([comp for seq_id, comp in consensus.comp_to_all_sequences.items() if seq_id in consensus.sequences]))
        if consensus.parent is not None:
            tree_graph.add_edge(consensus.parent, consensus.id, weight=len(consensus.sequences))
    return tree_graph

def get_consensus_tree_graph(jsonpangenome: JSONPangenome, tree, sliderValue):
    # read positions
    tree.graph.setdefault('graph', {})['rankdir'] = 'LR'
    node_id_to_x_y = graphviz_layout(tree, prog='dot')


    # prepare dots todo uprościć
    dots_labels = [tree.nodes[node_id] for node_id in range(len(node_id_to_x_y))]
    dots_labels_on_hover = [f'min_comp: {tree.nodes[node_id]["mincomp"]}' for node_id in range(len(node_id_to_x_y))]
    # dots_labels_on_hover = [f'min_comp: {tree.nodes[node_id]["mincomp"]}\nsequences: {tree.nodes[node_id]["sequences"]}' for node_id in range(len(node_id_to_x_y))]
    dots_numbers = [n for n in range(len(node_id_to_x_y))]
    dots_positions = [[tree.nodes[node_id]["mincomp"]*100, node_id_to_x_y[node_id][1]] for node_id in range(len(node_id_to_x_y))] #todo sprawdzic czy dobrze
    dots_x = [dot_x for [dot_x, _] in dots_positions] #todo czy to będzie dobrze posortowane?
    dots_y = [dot_y for [_, dot_y] in dots_positions] #todo czy to będzie dobrze posortowane?
    dots_annotations = [{'x':x_pos,
                         'y':y_pos,
                         'text':f"{i}",
                         'showarrow':False}
                        for i, (x_pos, y_pos) in enumerate(zip(dots_x, dots_y))]

    lines_x = []
    lines_y = []
    for u, v in tree.edges:
        lines_x += [dots_positions[u][0], dots_positions[v][0], None]
        lines_y += [dots_positions[u][1], dots_positions[v][1], None]
    # < class 'list'>:
    lines = go.Scatter(x=lines_x,
                       y=lines_y,
                       mode='lines',
                       line=dict(color='rgb(210,210,210)', width=1),
                       hoverinfo='none'
                       )
    line=go.Scatter(x=[sliderValue*100, sliderValue*100],
                    y=[0, 200],
                    mode='lines')
    dots = go.Scatter(x=dots_x,
                      y=dots_y,
                      mode='markers',
                      name='',
                      marker=dict(symbol='circle',
                                  size=25,
                                  color='white',
                                  line=dict(color='rgb(50,50,50)',
                                            width=1)
                                  ),
                      text=dots_labels_on_hover,
                      hoverinfo='text',
                      opacity=0.8,
                      customdata=dots_numbers#todo ?
                      )

    # graph settings - hide axis line, grid, ticklabels and  title
    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                )

    #graph settings - layout
    layout = dict(title='Consensuses Tree',
                  annotations=dots_annotations,
                  font=dict(size=12),
                  showlegend=False,
                  xaxis=go.layout.XAxis(axis),
                  yaxis=go.layout.YAxis(axis),
                  margin=dict(l=40, r=40, b=85, t=100),
                  hovermode='closest',
                  plot_bgcolor='rgb(248,248,248)',
                  width=1000,
                  height=700
                  )

    return go.Figure(
            data=[lines,dots, line],
            layout=layout
            )