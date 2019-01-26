import unittest
from ddt import ddt, data

# from context import mafreader
from context import metadatareader
from context import Node
from context import nucleotides as n
from context import Pangraph
from context import pathtools
from context import maf_to_dagmaf

from graph.Pangraph import PangraphBuilderFromDAG


@ddt
class MafreaderTest(unittest.TestCase):

    def setUp(self):
        self.test1_metadata = metadatareader.read(pathtools.get_file_content("Files/test1_metadata.json"))
        self.ebola_metadata = metadatareader.read(pathtools.get_file_content("Files/Ebola/ebola_metadata.json"))
        self.mycoplasma_metadata = metadatareader.read(pathtools.get_file_content("Files/Mycoplasma/mycoplasma_metadata.json"))

    @data("Files/test1.maf")
    # @unittest.skip("simple1")
    def test_read_maf1(self, maf_path):
        expected_nodes = [
            Node(id=0, base=n.code('C'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('T'), in_nodes=[], aligned_to=None),
            Node(id=2, base=n.code('G'), in_nodes=[1], aligned_to=None),
            Node(id=3, base=n.code('T'), in_nodes=[2], aligned_to=None),
            Node(id=4, base=n.code('G'), in_nodes=[3], aligned_to=5),
            Node(id=5, base=n.code('T'), in_nodes=[0], aligned_to=4),
            Node(id=6, base=n.code('A'), in_nodes=[4], aligned_to=None),
            Node(id=7, base=n.code('A'), in_nodes=[5], aligned_to=None),
            Node(id=8, base=n.code('C'), in_nodes=[6], aligned_to=None),
        ]

        expected_pats = {
            "testseq0": [1, 2, 3],
            "testseq1": [0, 5, 7],
            "testseq2": [3, 4, 6, 8]
        }
        expected_pangraph = Pangraph()
        expected_pangraph._nodes = expected_nodes
        expected_pangraph.set_paths(len(expected_nodes), expected_pats)
        dagmaf = maf_to_dagmaf(maf_path)
        actual_pangraph = Pangraph()
        PangraphBuilderFromDAG.build(dagmaf, actual_pangraph, self.test1_metadata)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    @data("Files/test2.maf")
    # @unittest.skip("simple2")
    def test_read_maf2(self, maf_path):
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
            Node(id=2, base=n.code('T'), in_nodes=[1], aligned_to=None),
            Node(id=3, base=n.code('G'), in_nodes=[1, 2], aligned_to=None),
            Node(id=4, base=n.code('G'), in_nodes=[3], aligned_to=None),
            Node(id=5, base=n.code('A'), in_nodes=[3], aligned_to=6),
            Node(id=6, base=n.code('T'), in_nodes=[4], aligned_to=5),
            Node(id=7, base=n.code('C'), in_nodes=[6], aligned_to=None)
        ]

        expected_pats = {
            "testseq1": [0, 1, 3, 4, 6, 7],
            "testseq2": [0, 1, 2, 3, 5],
        }

        expected_pangraph = Pangraph()
        expected_pangraph._nodes = expected_nodes
        expected_pangraph.set_paths(len(expected_nodes), expected_pats)
        dagmaf = maf_to_dagmaf(maf_path)
        actual_pangraph = Pangraph()
        PangraphBuilderFromDAG.build(dagmaf, actual_pangraph, self.test1_metadata)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


    @data("Files/test3.maf")
    @unittest.skip("empty")
    def test_read_maf3(self, maf_path):
        expected_nodes = []

        expected_pats = {
            "testseq0": [],
            "testseq1": []
        }

        expected_pangraph = Pangraph()
        expected_pangraph._nodes = expected_nodes
        expected_pangraph.set_paths(len(expected_nodes), expected_pats)
        dagmaf = maf_to_dagmaf(maf_path)
        actual_pangraph = Pangraph()
        PangraphBuilderFromDAG.build(dagmaf, actual_pangraph, self.test1_metadata)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    @data("Files/test4.maf")
    # @unittest.skip("simple4")
    def test_read_maf4(self, maf_path):
        expected_nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None)]

        expected_paths = {
            "testseq0": [0],
            "testseq1": [0],
            "testseq2": [0],
            "testseq3": [0]
        }

        expected_pangraph = Pangraph()
        expected_pangraph._nodes = expected_nodes
        expected_pangraph.set_paths(len(expected_nodes), expected_paths)
        dagmaf = maf_to_dagmaf(maf_path)
        actual_pangraph = Pangraph()
        PangraphBuilderFromDAG.build(dagmaf, actual_pangraph, self.test1_metadata)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


    @data("Files/test5.maf")
    # @unittest.skip("simple5")
    def test_read_maf5(self, maf_path):
        expected_nodes = [
                            Node(id=0, base=n.code('C'), in_nodes=[], aligned_to=None),
                            Node(id=1, base=n.code('T'), in_nodes=[], aligned_to=None),
                            Node(id=2, base=n.code('G'), in_nodes=[1], aligned_to=None),
                            Node(id=3, base=n.code('T'), in_nodes=[2], aligned_to=None),
                            Node(id=4, base=n.code('G'), in_nodes=[3], aligned_to=5),
                            Node(id=5, base=n.code('T'), in_nodes=[0], aligned_to=4),
                            Node(id=6, base=n.code('A'), in_nodes=[4], aligned_to=None),
                            Node(id=7, base=n.code('A'), in_nodes=[5], aligned_to=None),
                            Node(id=8, base=n.code('C'), in_nodes=[6], aligned_to=None)]

        expected_pats = {
            "testseq0": [1, 2, 3],
            "testseq1": [0, 5, 7],
            "testseq2": [3, 4, 6, 8]
        }

        expected_pangraph = Pangraph()
        expected_pangraph._nodes = expected_nodes
        expected_pangraph.set_paths(len(expected_nodes), expected_pats)
        dagmaf = maf_to_dagmaf(maf_path)
        actual_pangraph = Pangraph()
        PangraphBuilderFromDAG.build(dagmaf, actual_pangraph, self.test1_metadata)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


    @data("Files/test6.maf")
    def test_read_maf6(self, maf_path):
        expected_nodes = [
                            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
                            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
                            Node(id=2, base=n.code('T'), in_nodes=[0,1], aligned_to=None),
                            Node(id=3, base=n.code('G'), in_nodes=[], aligned_to=None),
                            Node(id=4, base=n.code('A'), in_nodes=[2], aligned_to=5),
                            Node(id=5, base=n.code('G'), in_nodes=[3], aligned_to=4),
                            Node(id=6, base=n.code('C'), in_nodes=[4], aligned_to=None),
                            Node(id=7, base=n.code('A'), in_nodes=[6], aligned_to=None),
                            Node(id=8, base=n.code('T'), in_nodes=[6,7], aligned_to=None),
                            Node(id=9, base=n.code('G'), in_nodes=[8], aligned_to=None),
                            Node(id=10, base=n.code('A'), in_nodes=[9], aligned_to=None),
                            Node(id=11, base=n.code('A'), in_nodes=[9,10], aligned_to=None),
                            Node(id=12, base=n.code('A'), in_nodes=[11], aligned_to=13),
                            Node(id=13, base=n.code('C'), in_nodes=[11], aligned_to=12),
                            Node(id=14, base=n.code('A'), in_nodes=[12], aligned_to=15),
                            Node(id=15, base=n.code('T'), in_nodes=[12,13], aligned_to=14),
        ]

        expected_pats = {
            "testseq1": [0, 1, 2, 4, 6, 8, 9, 10, 11, 12, 14],
            "testseq2": [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 13, 15],
            "testseq3": [0, 2, 3, 5, 6, 7, 8, 9, 11, 12, 15]
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
        PangraphBuilderFromDAG.build(dagmaf, pangraph, self.test1_metadata)
        return pangraph

    @data("Files/Ebola/ebola.maf")
    @unittest.skip("ebola")
    def test_read_maf_ebola(self, maf_path):
        pangraph = mafreader.read(maf_path, self.ebola_metadata)

    @data("Files/Mycoplasma/mycoplasma.maf")
    @unittest.skip("myco")
    def test_read_maf_myco(self, maf_path):
        pangraph = mafreader.read(maf_path, self.mycoplasma_metadata)

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
