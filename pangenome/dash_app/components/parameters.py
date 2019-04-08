from typing import Dict

from pangenome.pang.fileformats.json.JSONPangenome import JSONPangenome


def get_data(jsonpangenome: JSONPangenome) -> Dict[str, str]:
    parameters_data = {}
    if jsonpangenome['program_parameters']:
        parameters_data["Multialignment file name"] = jsonpangenome.program_parameters.multialignment_file_path
        parameters_data["Datatype"] = jsonpangenome.program_parameters.datatype
        parameters_data["Metadata file name"] = jsonpangenome.program_parameters.metadata_file_path
        parameters_data["Blosum file name"] = jsonpangenome.program_parameters.blosum_file_path

        parameters_data["Output po"] = "Yes" if jsonpangenome.program_parameters.output_po else "No"
        parameters_data["Generate fasta"] = "Yes" if jsonpangenome.program_parameters.generate_fasta else "No"
        parameters_data["Include nodes in output"] = "Yes" if jsonpangenome.program_parameters.output_with_nodes else "No"

        parameters_data["Build from DAG"] = "No" if jsonpangenome.program_parameters.raw_maf else "Yes"
        parameters_data["Fasta complementation option"] = jsonpangenome.program_parameters.fasta_complementation_option
        parameters_data["Missing base symbol"] = jsonpangenome.program_parameters.missing_base_symbol
        parameters_data["Fasta source file"] = jsonpangenome.program_parameters.fasta_source_file

        parameters_data["Consensus type"] = str(jsonpangenome.program_parameters.consensus_type)
        parameters_data["HBMIN"] = float(jsonpangenome.program_parameters.hbmin)

        parameters_data["MAX Cutoff Strategy"] = str(jsonpangenome.program_parameters.max_cutoff_strategy)
        parameters_data["RANGE"] = str(jsonpangenome.program_parameters.search_range)
        parameters_data["NODE Cutoff Strategy"] = str(jsonpangenome.program_parameters.node_cutoff_strategy)
        parameters_data["MULTIPLIER"] = float(jsonpangenome.program_parameters.multiplier)
        parameters_data["STOP"] = float(jsonpangenome.program_parameters.stop)
        parameters_data["RE CONSENSUS"] = jsonpangenome.program_parameters.re_consensus
        parameters_data["P"] = jsonpangenome.program_parameters.p

    parameters_data["Nodes count"] = len(jsonpangenome.nodes) if jsonpangenome.nodes else 0
    parameters_data["Sequences count"] = len(jsonpangenome.sequences) if jsonpangenome.consensuses else 0
    if jsonpangenome.consensuses:
        parameters_data["Consensus tree"] = f"Consensus tree with {len(jsonpangenome.consensuses)} nodes generated."
    else:
        parameters_data["Consensus tree"] = f"No consensus tree generated."
    return parameters_data


