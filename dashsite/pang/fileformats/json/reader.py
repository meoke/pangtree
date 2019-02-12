import jsonpickle


def json_to_jsonpangenome(pangenome_json):
    return jsonpickle.decode(pangenome_json)
