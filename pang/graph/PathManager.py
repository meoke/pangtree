import numpy as np
from typing import Dict, List

PathsNodesDict = Dict[str, int]
PathsNamesList = List[str]
#todo czy path to aby na pewno string? a nie path/sequence?

class PathManager:
    def __init__(self, start_node_id: int=0, max_nodes_count: int=0, paths_names: PathsNamesList=None):
        # if paths_to_node_ids:
        #     self.init_from_dict(paths_to_node_ids)
        # else:
        paths_names = [] if not paths_names else paths_names
        self.paths = np.zeros(shape=(len(paths_names), max_nodes_count), dtype=bool)
        self.path_names_to_array_id = {path_name: i for i, path_name in enumerate(sorted(paths_names))}
        self.start_node_id = start_node_id

    def init_from_dict(self, paths_to_node_ids):
        max_node_id = max([node_id for node_ids in paths_to_node_ids.values() for node_id in node_ids] )
        self.start_node_id = 0
        self.paths = np.zeros(shape=(len(paths_to_node_ids.keys()), max_node_id+1), dtype=bool)
        self.path_names_to_array_id = {path_name: array_id
                                       for array_id, path_name in enumerate(sorted(paths_to_node_ids.keys()))}
        for path, node_ids in paths_to_node_ids.items():
            for node_id in node_ids:
                self.paths[self.path_names_to_array_id[path], node_id] = True

    def __eq__(self, other):
        return (np.array_equal(self.paths, other.paths)
                and self.path_names_to_array_id == other.path_names_to_array_id
                and self.start_node_id == other.start_node_id)
    #
    # def no_paths(self):
    #     return not any(x for p in self.paths for x in p)

    def mark(self, path_name, node_id):
        array_id = self.path_names_to_array_id[path_name]
        self.paths[array_id, node_id-self.start_node_id] = True

    def update(self, pathmanager, start):
        for path_name, array_id in pathmanager.path_names_to_array_id.items():
            current_array_id = self.path_names_to_array_id[path_name]
            end = start+np.shape(pathmanager.paths)[1]
            self.paths[current_array_id, start:end] = pathmanager.paths[array_id, :]

    def trim(self, nodes_count):
        self.paths = np.delete(self.paths, np.s_[nodes_count:], 1)

    def get_in_nodes(self, node_id):
        in_nodes = []
        node_ids_to_check = range(node_id-1, -1, -1)
        for seq in self.paths:
            if seq[node_id]:
                for col_id in node_ids_to_check:
                    if seq[col_id]:
                        in_nodes.append(col_id)
                        break
        return list(set(in_nodes))

    def get_paths_count(self):
        return self.paths.shape[0]

    def get_start_node_id(self, pathname):
        return self.get_nodes(pathname)[0]

    def get_sources_weights(self):
        unweighted_sources_weights = []
        for path in self.path_names_to_array_id.keys():
            path_nodes_ids = self.get_nodes(path)
            average_paths_in_nodes_count = self.get_average_paths_in_nodes_count(path_nodes_ids)
            unweighted_sources_weights.append(average_paths_in_nodes_count)

        max_weight = max(unweighted_sources_weights)
        min_weight = min(unweighted_sources_weights)
        diff_weight = max_weight - min_weight
        if diff_weight == 0:
            normalized_sources_weights = [100 for _ in range(len(unweighted_sources_weights))]
        else:
            normalized_sources_weights = [int((weight - min_weight)/diff_weight*100) for weight in unweighted_sources_weights]
        return normalized_sources_weights

    def get_sources_ids(self, node_id: int) -> List[int]:
        return (np.where(self.paths[:, node_id])[0]).tolist()

    def get_average_paths_in_nodes_count(self, path_nodes_ids):
        return np.mean(sum(self.paths[:, path_nodes_ids]))

    def get_nodes(self, pathname):
        source_id = self.path_names_to_array_id[pathname]
        return np.where(self.paths[source_id, :])[0]

    def get_path_nodes_count(self, pathname):
        return len(self.get_nodes(pathname))

