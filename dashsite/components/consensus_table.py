from collections import deque
import numpy as np

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
    consensuses_names =  [consensus_tree.nodes[node_id]['name'] for node_id in consensuses_nodes_ids]
    data=np.array([['SequenceID', 'Name', 'Title', 'Group'] + consensuses_names])
    for seq in jsonpangenome.sequences:
        row = np.array([[seq.id, seq.name, 'title', 'grupa'] + \
              [consensus_tree.nodes[node_id]['comp'][seq.name] for node_id in consensuses_nodes_ids]])
        data = np.append(data, row, axis=0)
    return data


def get_table_rows(jsonified_consensuses_table_data):
    def convert_to_number(s):
        try:
            return float(s)
        except ValueError:
            try:
                return int(s)
            except:
                return s

    first_item = jsonified_consensuses_table_data[0]

    rows = []
    for row in jsonified_consensuses_table_data[1:]:
        dict_row = {v: convert_to_number(row[k]) for k, v in first_item.items()}
        rows.append(dict_row)
    return rows


def get_table_columns(jsonified_consensuses_table_data):
    column_headers = jsonified_consensuses_table_data.to_dict().keys()
    return [{"name": c, "id": c} for i, c in enumerate(column_headers)]


