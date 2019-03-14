from typing import Dict

from fileformats.json.JSONPangenome import JSONPangenome


def get_data(jsonpangenome: JSONPangenome) -> Dict[str, str]:
    parameters_data = {}
    if jsonpangenome.program_parameters:
        parameters_data["Build from DAG"] = "No" if jsonpangenome.program_parameters.not_dag else "Yes"
        parameters_data["Fasta complementation option"] = str(jsonpangenome.program_parameters.fasta_complementation_option)

        parameters_data["Consensus type"] = str(jsonpangenome.program_parameters.consensus_type)
        parameters_data["MAX Cutoff Strategy"] = str(jsonpangenome.program_parameters.max_cutoff_strategy)
        parameters_data["NODE Cutoff Strategy"] = str(jsonpangenome.program_parameters.node_cutoff_strategy)
        parameters_data["HBMIN"] = float(jsonpangenome.program_parameters.hbmin)
        parameters_data["RANGE"] = str(jsonpangenome.program_parameters.r)
        parameters_data["MULTIPLIER"] = float(jsonpangenome.program_parameters.multiplier)
        parameters_data["STOP"] = float(jsonpangenome.program_parameters.stop)
        parameters_data["RE CONSENSUS"] = jsonpangenome.program_parameters.re_consensus

    parameters_data["Nodes count"] = len(jsonpangenome.nodes) if jsonpangenome.nodes else 0
    parameters_data["Sequences count"] = len(jsonpangenome.sequences) if jsonpangenome.consensuses else 0
    if jsonpangenome.consensuses:
        parameters_data["Consensus tree"] = f"Consensus tree with {len(jsonpangenome.consensuses)} nodes generated."
    else:
        parameters_data["Consensus tree"] = f"No consensus tree generated."
    return parameters_data


