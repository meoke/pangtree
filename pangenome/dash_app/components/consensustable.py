from collections import deque
from typing import List, Dict

from consensus.ConsensusNode import ConsensusNodeID
from dash_app.components import consensustree
from dash_app.layout.css_styles import colors
from fileformats.json.JSONPangenome import JSONPangenome
import pandas as pd
import networkx as nx

def get_full_table_data(jsonpangenome: JSONPangenome) -> pd.DataFrame:
    if not jsonpangenome.sequences:
        return pd.DataFrame()
    if jsonpangenome.consensuses:
        all_consensuses_ids = [get_consensus_column_name(c.node_id) for c in jsonpangenome.consensuses]
    else:
        all_consensuses_ids = []
    first_consensus = jsonpangenome.sequences[0]
    consensus_df = pd.DataFrame(columns=["ID", "SEQID"] +
                                   [*first_consensus.metadata.keys()] +
                                   all_consensuses_ids)
    for seq in jsonpangenome.sequences:
        row = {"ID": seq.sequence_int_id,
               "SEQID": seq.sequence_str_id}
        for m, v in seq.metadata.items():
            row[m] = v
        for c in jsonpangenome.consensuses:
            row[get_consensus_column_name(c.node_id)] = get_mapped_compatibility(c.comp_to_all_sequences[seq.sequence_str_id])
        consensus_df = consensus_df.append(row, ignore_index=True)
    return consensus_df

def get_mapped_compatibility(compatibility: float) -> str:
    return "{:.{}f}".format(compatibility, 4)

def get_consensus_column_name(consensus_id: ConsensusNodeID) -> str:
    return f"CONSENSUS{consensus_id}"

def get_metadata_list(full_consensustable: pd.DataFrame) -> List[str]:
    return [colname for colname in list(full_consensustable) if "CONSENSUS" not in colname]

def remove_smaller_than_slider(full_consensustable_data: pd.DataFrame, tree: nx.DiGraph, slider_value: float) -> pd.DataFrame:
    consensuses_ids_to_hide = []
    nodes_to_visit = deque([0])
    while nodes_to_visit:
        current_node_id = nodes_to_visit.pop()
        current_node = tree.nodes[current_node_id]

        if current_node['mincomp'] > slider_value:
            current_node_offspring_nodes_ids = consensustree.get_offspring_ids(tree, current_node_id)
            consensuses_ids_to_hide += current_node_offspring_nodes_ids
        else:
            consensuses_ids_to_hide.append(current_node_id)
            nodes_to_visit.extend(current_node['children_consensuses'])
    columns_to_hide = map(lambda c_id : get_consensus_column_name(c_id), consensuses_ids_to_hide)

    return full_consensustable_data.drop(columns_to_hide, axis=1)


def hide_children(consensus_tree, parent_id):
    nodes_to_visit = deque(consensus_tree.nodes[parent_id]['children_consensuses'])
    while nodes_to_visit:
        current_node_id = nodes_to_visit.pop()
        consensus_tree.nodes[current_node_id]['show_in_table'] = False
        nodes_to_visit.extend(consensus_tree.nodes[current_node_id]['children_consensuses'])
    return consensus_tree


def get_consensus_details_df(consensus_node_id: ConsensusNodeID,
                             full_consensustable: pd.DataFrame,
                             tree: nx.DiGraph) -> pd.DataFrame:
    sequences_ids = tree.nodes[consensus_node_id]['sequences_ids']
    columns_names_to_show = get_metadata_list(full_consensustable) + [get_consensus_column_name(consensus_node_id)]
    return full_consensustable[columns_names_to_show].loc[full_consensustable["ID"].isin(sequences_ids)]


def get_cells_styling(tree: nx.DiGraph, partial_consensustable_data: pd.DataFrame) -> List[Dict]:

    consensuses_colnames = [colname for colname in list(partial_consensustable_data) if "CONSENSUS" in colname]
    styling_conditions = []
    for consensus_colname in consensuses_colnames:
        consensus_id = int(consensus_colname[9:])
        consensus_mincomp = tree.nodes[consensus_id]['mincomp']
        styling_conditions.append(get_cell_styling_dict(consensus_colname, get_mapped_compatibility(consensus_mincomp)))
    return styling_conditions

def get_cell_styling_dict(consensus_colname, mincomp):
    return {
        'if': {'column_id': f'{consensus_colname}',
               'filter': f'{consensus_colname} >= "{mincomp}"'},
        'backgroundColor': colors['warm_background']}
