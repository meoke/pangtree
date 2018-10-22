import numpy as np
from typing import Dict, List
from graph.errors import NoPath

PathsNodesDict = Dict[str, int]
PathsNamesList = List[str]
#todo czy path to aby na pewno string? a nie path/sequence?


class PathManager:
    def __init__(self, start_node_id: int = 0, max_nodes_count: object = 0, paths_names: object = None) -> object:
        paths_names = [] if not paths_names else paths_names
        self.paths = np.zeros(shape=(len(paths_names), max_nodes_count), dtype=bool)
        self.path_names_to_array_id = {path_name: i for i, path_name in enumerate(sorted(paths_names))}
        self.start_node_id = start_node_id

    def init_from_dict(self, paths_to_node_ids):
        if not paths_to_node_ids.keys():
            return
        max_node_id = max([node_id for node_ids in paths_to_node_ids.values() for node_id in node_ids])
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

    def get_sources_weights_dict(self):
        unweighted_sources_weights = {}
        for path in self.path_names_to_array_id.keys():
            path_nodes_ids = self.get_nodes(path)
            average_paths_in_nodes_count = self.get_average_paths_in_nodes_count(path_nodes_ids)
            unweighted_sources_weights[path] = average_paths_in_nodes_count

        max_weight = max(unweighted_sources_weights.values())
        min_weight = min(unweighted_sources_weights.values())
        diff_weight = max_weight - min_weight
        if diff_weight == 0:
            normalized_sources_weights_dict = {path_key: 100 for path_key in unweighted_sources_weights.keys()}
        else:
            normalized_sources_weights_dict = {path: int((weight - min_weight)/diff_weight*100) for path, weight in unweighted_sources_weights.items()}
        return normalized_sources_weights_dict

    def get_sources_ids(self, node_id: int) -> List[int]:
        return (np.where(self.paths[:, node_id])[0]).tolist()

    def get_average_paths_in_nodes_count(self, path_nodes_ids):
        if len(path_nodes_ids):
            return np.mean(sum(self.paths[:, path_nodes_ids]))
        return 0

    def get_nodes(self, pathname):
        source_id = self.path_names_to_array_id[pathname]
        return np.where(self.paths[source_id, :])[0]

    def get_path_nodes_count(self, pathname):
        return len(self.get_nodes(pathname))

    def get_path_names(self):
        return list(self.path_names_to_array_id.keys())

    def get_path_ids(self):
        return list(self.path_names_to_array_id.values())

    def keep_paths_ids(self, sequences_ids: List[int]):
        """Removes unnecessary paths info and return ids of still necessary nodes."""
        # paths_ids_to_delete = list(set(self.path_names_to_array_id.values()) - set(sequences_ids))
        # self.paths = np.delete(self.paths, paths_ids_to_delete, 0)
        # return np.where(np.logical_or.reduce(self.paths))[0]
        path_names_to_delete = []
        path_ids_to_delete = []
        d_copy = self.path_names_to_array_id.copy()
        i = 0
        for path_name, path_id in d_copy.items():
            if path_id in sequences_ids:
                self.path_names_to_array_id[path_name] = i
                i += 1
            else:
                path_names_to_delete.append(path_name)
                path_ids_to_delete.append(path_id)
        self.paths = np.delete(self.paths, path_ids_to_delete, 0)
        for path_name in path_names_to_delete:
            del self.path_names_to_array_id[path_name]

    def keep_paths_names(self, sequences_names: List[str]):
        """Removes unnecessary paths info and return ids of still necessary nodes."""
        # paths_ids_to_delete = list(set(self.path_names_to_array_id.values()) - set(sequences_ids))
        # self.paths = np.delete(self.paths, paths_ids_to_delete, 0)
        # return np.where(np.logical_or.reduce(self.paths))[0]
        path_names_to_delete = []
        path_ids_to_delete = []
        d_copy = self.path_names_to_array_id.copy()
        i = 0
        for path_name, path_id in d_copy.items():
            if path_name in sequences_names:
                self.path_names_to_array_id[path_name] = i
                i += 1
            else:
                path_names_to_delete.append(path_name)
                path_ids_to_delete.append(path_id)
        self.paths = np.delete(self.paths, path_ids_to_delete, 0)
        for path_name in path_names_to_delete:
            del self.path_names_to_array_id[path_name]

    def keep_nodes_ids(self, nodes_ids_to_keep):
        inactive_nodes = [node_id for node_id in range(self.paths.shape[1]) if node_id not in nodes_ids_to_keep]
        self.paths = np.delete(self.paths, inactive_nodes, 1)

    def get_top_consensus(self):
        return self.paths[0]

    def remove_empty_paths(self):
        paths_names_to_keep = [path_id
                               for path_name, path_id
                               in self.path_names_to_array_id.items()
                               if any(self.paths[path_id])]
        self.keep_paths_ids(paths_names_to_keep)

    def get_active_nodes(self):
        return np.where(np.logical_or.reduce(self.paths))[0]

    def get_paths(self):
        return self.paths

    def get_path_name(self, path_id):
        pathname = [pathname for pathname, array_id in self.path_names_to_array_id.items() if array_id == path_id]
        if len(pathname) == 0:
            raise NoPath
        return pathname[0]

    def add_path(self, pathname_prefix, path):
        path_id = self.paths.shape[0]
        self.path_names_to_array_id[f"{pathname_prefix}_{path_id}"] = path_id
        self.paths = np.append(self.paths, [path], axis=0)
        return path_id

    def get_path_id(self, pathname):
        return self.path_names_to_array_id[pathname]

    def clear_paths(self):
        self.path_names_to_array_id = {}
        self.paths = np.delete(self.paths, slice(0, self.paths.shape[0]), axis=0)

