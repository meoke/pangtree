import numpy as np


class PathManager:
    def __init__(self, start_node_id, max_nodes_count, paths_names):
        self.paths = np.zeros(shape=(len(paths_names), max_nodes_count), dtype=bool)
        self.path_names_to_array_id = {path_name: i for i, path_name in enumerate(paths_names)}
        self.start_node_id = start_node_id

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
