import unittest
from ddt import ddt, data

# from context import mafreader
from context import metadatareader
from context import Node
from context import nucleotides as n
from context import Pangraph
from context import pathtools
from context import maf_to_dagmaf
from context import FastaSource

from graph.Pangraph import PangraphBuilderFromDAG


@ddt
class FillMafGapsTest(unittest.TestCase):
    class FakeFastaSource(FastaSource):
        def __init__(self):
            self.testseq0 = ""
            self.testseq1 = "ACTAGGT"
            self.testseq2 = "GGTCAGT"
            self.testseq3 = ""
            self.testseq4 = ""

        def get_source(self, id: str, start: int = None, end: int = None) -> str:
            if id == "testseq0":
                return self.testseq0[start: end]
            elif id == "testseq1":
                return self.testseq1[start: end]
            elif id == "testseq2":
                return self.testseq2[start: end]
            elif id == "testseq3":
                return self.testseq3[start: end]
            elif id == "testseq4":
                return self.testseq4[start: end]
            else:
                raise Exception("No record found with given id!")

    def setUp(self):
        self.test_metadata = metadatareader.read(pathtools.get_file_content("Files/test1_metadata.json"))
        self.fasta_source = FillMafGapsTest.FakeFastaSource()

    @data("Files/maf_gaps/test_1_left.maf")
    def test_maf_gap_1_left(self, maf_path):
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
            Node(id=2, base=n.code('T'), in_nodes=[1], aligned_to=None),
            Node(id=3, base=n.code('A'), in_nodes=[2], aligned_to=4),
            Node(id=4, base=n.code('G'), in_nodes=[], aligned_to=3),
            Node(id=5, base=n.code('G'), in_nodes=[3, 4], aligned_to=None),
            Node(id=6, base=n.code('G'), in_nodes=[5], aligned_to=7),
            Node(id=7, base=n.code('T'), in_nodes=[5], aligned_to=6),
            Node(id=8, base=n.code('C'), in_nodes=[7], aligned_to=None),
            Node(id=9, base=n.code('A'), in_nodes=[8], aligned_to=None),
            Node(id=10, base=n.code('G'), in_nodes=[9], aligned_to=None),
            Node(id=11, base=n.code('T'), in_nodes=[6, 10], aligned_to=None),
        ]

        expected_pats = {
            "testseq1": [0, 1, 2, 3, 5, 6, 11],
            "testseq2": [4, 5, 7, 8, 9, 10, 11]
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
        builder = PangraphBuilderFromDAG(self.test_metadata, FillMafGapsTest.FakeFastaSource)
        builder.build(dagmaf, pangraph)
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
