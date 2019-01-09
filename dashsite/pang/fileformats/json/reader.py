import jsonpickle


def json_to_jsonpangenome(pangenome_json):
    return jsonpickle.decode(pangenome_json)


# todo to raczej niepotrzebne
# def json_to_pangenome(pangenome_json)-> Pangenome:
#     #todo type pangenome_json
#     jsonpangenome = json_to_jsonpangenome(pangenome_json)
#     return jsonpangenome.convert_to_pangenome()