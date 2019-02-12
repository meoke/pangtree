import unittest

from .context import Pangraph


class PangraphComparator(unittest.TestCase):
    def compare_pangraphs(self, actual_pangraph, expected_pangraph):
        try:
            self.assertEqual(actual_pangraph, expected_pangraph)
        except Exception as ex:
            self.show_pangraph_differences(actual_pangraph, expected_pangraph)
            raise ex

    def setup_pangraph(self, expected_nodes, expected_paths):
        pangraph = Pangraph()
        pangraph._nodes = expected_nodes
        pangraph.set_paths(len(expected_nodes), expected_paths)
        return pangraph

    def show_pangraph_differences(self, actual_pangraph, expected_pangraph):
        print("Nodes differences: ")
        self.show_nodes_differences(actual_pangraph._nodes,
                                                  expected_pangraph._nodes)

        print("Paths differences: ")
        self.show_paths_differences(actual_pangraph._pathmanager.paths,
                                                  expected_pangraph._pathmanager.paths)

    def show_nodes_differences(self, actual_nodes, expected_nodes):
        if len(actual_nodes) != len(expected_nodes):
            print(f"Actual graph has {len(actual_nodes)} nodes while expected graph: {len(expected_nodes)}")
        for node_id in range(len(actual_nodes)):
            if actual_nodes[node_id] is None or actual_nodes[node_id] != expected_nodes[node_id]:
                print(f"Nodes {node_id} differ:")
                print(f"Actual: {actual_nodes[node_id]}")
                print(f"Expected: {expected_nodes[node_id]}")

    def show_paths_differences(self, actual_paths, expected_paths):
        if actual_paths.shape != expected_paths.shape:
            print(f"Actual pm has shape {actual_paths.shape} while expected pm: {expected_paths.shape}")

        for i, (actual_row, expected_row) in enumerate(zip(actual_paths, expected_paths)):
            if (actual_row != expected_row).any():
                print(f"Rows {i} differ.")
                print(f"Actual row: {actual_row}")
                print(f"Expected row: {expected_row}")
