import pandas as pd


def get_details(tree_click_data, tree, jsonpangenome):
    clicked_node = tree_click_data['points'][0]
    node_id = clicked_node['pointIndex']
    sequences_ids = tree.nodes[node_id]['sequences_ids']

    sequences_data = pd.DataFrame(columns=['ID', 'Name', 'Group'])
    for seq_id in sequences_ids:
        row = {"ID": seq_id,
               "Name": jsonpangenome.sequences[seq_id].name,
               "Group": jsonpangenome.sequences[seq_id].group}
        sequences_data = sequences_data.append(row, ignore_index=True)
    return sequences_data.to_dict("rows")
