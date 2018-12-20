from collections import deque
import pandas as pd


def mark_nodes_to_show(consensus_tree, slider_value):
    nodes_to_visit = deque([0])
    while nodes_to_visit:
        current_node_id = nodes_to_visit.pop()
        current_node = consensus_tree.nodes[current_node_id]
        if current_node['mincomp'] > slider_value:
            current_node['show_in_table'] = True
            consensus_tree = hide_children(consensus_tree, current_node_id)
        else:
            current_node['show_in_table'] = False
            nodes_to_visit.extend(current_node['children_consensuses'])
    return consensus_tree


def hide_children(consensus_tree, parent_id):
    nodes_to_visit = deque(consensus_tree.nodes[parent_id]['children_consensuses'])
    while nodes_to_visit:
        current_node_id = nodes_to_visit.pop()
        consensus_tree.nodes[current_node_id]['show_in_table'] = False
        nodes_to_visit.extend(consensus_tree.nodes[current_node_id]['children_consensuses'])
    return consensus_tree


def get_consensus_table_data(jsonpangenome, consensus_tree, slider_value):
    consensus_tree = mark_nodes_to_show(consensus_tree, slider_value)
    consensuses_nodes_ids = [n for n in consensus_tree.nodes if consensus_tree.nodes[n]['show_in_table']]
    consensuses_names = [consensus_tree.nodes[node_id]['name'] for node_id in consensuses_nodes_ids]
    df = pd.DataFrame(columns=['SequenceID', 'Name', 'Title', 'Group'] + consensuses_names)
    for seq in jsonpangenome.sequences:
        row = {}
        for node_id in consensuses_nodes_ids:
            consensus_name = consensus_tree.nodes[node_id]['name']
            row[consensus_name] = float(consensus_tree.nodes[node_id]['comp'][seq.name])
        row["SequenceID"] = int(seq.id)
        row["Name"] = seq.name
        row["Title"] = 'title'
        row["Group"] = seq.group
        df = df.append(row, ignore_index=True)
    return df


def s(df_table_data, tree):
    df_json = df_table_data.to_dict("rows")
    # df_json["10"]["C_0"]
    return df_json



