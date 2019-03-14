from typing import List

from consensus.ConsensusNode import ConsensusNodeID
from fileformats.json.JSONPangenome import JSONPangenome
import pandas as pd


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
            row[get_consensus_column_name(c.node_id)] = c.comp_to_all_sequences[seq.sequence_str_id]
        consensus_df = consensus_df.append(row, ignore_index=True)
    for consensus_id in all_consensuses_ids:
        consensus_df[consensus_id] = consensus_df[consensus_id].map('{:,.3f}'.format).map(str)

    return consensus_df

def get_consensus_column_name(consensus_id: ConsensusNodeID) -> str:
    return f"CONSENSUS{consensus_id}"

def get_metadata_list(full_consensustable: pd.DataFrame) -> List[str]:
    return [colname for colname in list(full_consensustable) if "CONSENSUS" not in colname]