from collections import deque
import pandas as pd
from app_style import colors


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



def table_to_rows_json(df_table_data):
    df_json = df_table_data.to_dict("rows")
    return df_json


def get_cells_styling(consensus_tree, table_data):
    consensuses_names = [header for header in list(table_data) if header[0:9] == 'CONSENSUS']
    styling_conditions = []
    for c in consensuses_names:
        consensus_id = int(c[9:])
        consensus_mincomp = consensus_tree.nodes[consensus_id]['mincomp']
        styling_conditions.append(get_cell_styling_dict(c, consensus_mincomp))
    return styling_conditions


def get_cell_styling_dict(consensus_name, mincomp):
    return {
        'if': {'column_id': f'{consensus_name}', 'filter': f'{consensus_name} >= "{mincomp}"'},
        'backgroundColor': colors['warm_background']}


def get_consensus_table_data(jsonpangenome, consensus_tree, slider_value) -> pd.DataFrame:
    consensus_tree = mark_nodes_to_show(consensus_tree, slider_value)
    return get_consensuses_table(jsonpangenome, consensus_tree)


def get_consensuses_table(jsonpangenome, consensus_tree):
    consensuses_nodes_ids = [n for n in consensus_tree.nodes if consensus_tree.nodes[n]['show_in_table']]
    consensuses_names = [consensus_tree.nodes[node_id]['name'] for node_id in consensuses_nodes_ids]
    table_data = pd.DataFrame(
        columns=['ID', 'Genbank ID', 'Assembly ID', 'Mafname', 'Name', 'Group', 'Leaf assignment'] +
                consensuses_names)

    for seq in jsonpangenome.sequences:
        leaf = [node_id
                for node_id
                in consensus_tree.nodes
                if seq.id in consensus_tree.nodes[node_id]['sequences_ids']
                and not consensus_tree.nodes[node_id]['children_consensuses']][0]
        row = {"ID": seq.id,
               "Genbank ID": seq.genbankID,
               "Assembly ID": seq.assemblyID,
               "Mafname": seq.mafname,
               "Name": seq.name,
               "Group": seq.group,
               "Leaf assignment": leaf}
        for node_id in consensuses_nodes_ids:
            consensus_name = consensus_tree.nodes[node_id]['name']
            row[consensus_name] = consensus_tree.nodes[node_id]['comp'][seq.mafname]
        table_data = table_data.append(row, ignore_index=True)
    for consensus_name in consensuses_names:
        table_data[consensus_name] = table_data[consensus_name].map('{:,.4f}'.format).map(str)
    return table_data


def get_all_consensuses_table(jsonpangenome, consensuses_tree):
    for node_id in consensuses_tree.nodes:
        consensuses_tree.nodes[node_id]['show_in_table'] = True
    return get_consensuses_table(jsonpangenome, consensuses_tree)