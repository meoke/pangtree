import unittest
from ddt import ddt, data

# from context import mafreader
from context import metadatareader
from context import Node
from context import nucleotides as n
from context import Pangraph
from context import pathtools
from context import maf_to_dagmaf
from context import FastaFileSystemSource

from tests_fillmafgaps import FillMafGapsTest

from graph.Pangraph import PangraphBuilderFromDAG


@ddt
class PangraphBuilderFastaSourceTest(unittest.TestCase):

    def setUp(self):
        self.fasta_source_dir = "Files/seq/fasta"
        self.seq_metadata = metadatareader.read(pathtools.get_file_content("Files/seq/seq_metadata.json"))

    @data("Files/seq/test9_parallel_first_blocks.maf")
    def test_read_maf0(self, maf_path):
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
            Node(id=2, base=n.code('T'), in_nodes=[1], aligned_to=None),
            Node(id=3, base=n.code('G'), in_nodes=[2], aligned_to=None),
            Node(id=4, base=n.code('A'), in_nodes=[3], aligned_to=None),
            Node(id=5, base=n.code('C'), in_nodes=[4], aligned_to=None),
            Node(id=6, base=n.code('T'), in_nodes=[5], aligned_to=None),
            Node(id=7, base=n.code('G'), in_nodes=[6], aligned_to=None),
            Node(id=8, base=n.code('A'), in_nodes=[7], aligned_to=None),
            Node(id=9, base=n.code('A'), in_nodes=[8], aligned_to=None),
        ]

        expected_pats = {
            "testseq0": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            "testseq1": [0, 1, 2, 3, 4, 5, 6, 7, 8],
        }
        expected_pangraph = self.setup_pangraph(expected_nodes, expected_pats)
        actual_pangraph = self.setup_pangraph_from_maf(maf_path)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def setup_pangraph(self, expected_nodes, expected_paths):
        pangraph = Pangraph()
        pangraph._nodes = expected_nodes
        pangraph.set_paths(len(expected_nodes), expected_paths)
        return pangraph

    def setup_pangraph_from_maf(self, maf_path):
        dagmaf = maf_to_dagmaf(maf_path)
        pangraph = Pangraph()
        builder_from_maf = PangraphBuilderFromDAG(self.test1_metadata, FastaFileSystemSource(self.fasta_source_dir))
        builder_from_maf.build(dagmaf, pangraph)
        return pangraph

    def compare_graphs(self, actual_graph, expected_graph):
        try:
            self.assertEqual(actual_graph, expected_graph)
        except Exception as ex:
            self.show_graph_differences(actual_graph=actual_graph, expected_graph=expected_graph)
            raise ex

    def compare_pangraphs(self, actual_pangraph, expected_pangraph):
        try:
            self.assertEqual(actual_pangraph, expected_pangraph)
        except Exception as ex:
            self.show_pangraph_differences(actual_pangraph, expected_pangraph)
            raise ex

    def show_pangraph_differences(self, actual_pangraph, expected_pangraph):
        print("Nodes differences: ")
        self.show_nodes_differences(actual_pangraph._nodes, expected_pangraph._nodes)

        print("Paths differences: ")
        self.show_paths_differences(actual_pangraph._pathmanager.paths,
                               expected_pangraph._pathmanager.paths)

    def show_nodes_differences(self, actual_nodes, expected_nodes):
        if len(actual_nodes) != len(expected_nodes):
            print(f"Actual graph has {len(actual_nodes)} nodes while expected graph: {len(expected_nodes)}")
        for node_id in range(len(actual_nodes)):
            if actual_nodes[node_id] != expected_nodes[node_id]:
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


if __name__ == '__main__':
    unittest.main()
