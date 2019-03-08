import unittest
from io import StringIO
from pathlib import Path

from tests.context import maf_to_dagmaf
from tests.context import Pangraph
from tests.context import PangraphBuilderFromDAG
from Bio import AlignIO

from tests.context import PangraphBuilderFromMAF


class PangraphBuilderTests(unittest.TestCase):
    @staticmethod
    def setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, metadata, fasta_source):
        mafcontent = PangraphBuilderTests.get_file_content(maf_path)
        pangraph = Pangraph()
        builder = PangraphBuilderFromDAG(metadata, fasta_source)
        builder.build(mafcontent, pangraph)
        return pangraph

    @staticmethod
    def setup_pangraph_from_maf(maf_path, metadata):
        mafalignment = PangraphBuilderTests.get_file_content(maf_path)
        pangraph = Pangraph()
        builder = PangraphBuilderFromMAF(metadata)
        builder.build(mafalignment, pangraph)
        return pangraph

    @staticmethod
    def setup_pangraph(expected_nodes, expected_paths):
        pangraph = Pangraph()
        pangraph.nodes = expected_nodes
        pangraph.paths = expected_paths
        return pangraph

    @staticmethod
    def get_file_content(path: Path) -> StringIO:
        """Returns file content."""

        with open(path) as input_file:
            return StringIO(input_file.read())

    @staticmethod
    def read_maf(path: Path):
        """Returns file content as BioPython AlignIO object."""

        return [*AlignIO.parse(path, "maf")]

    def compare_pangraphs(self, actual_pangraph, expected_pangraph):
        try:
            self.assertEqual(actual_pangraph, expected_pangraph)
        except Exception as ex:
            self.show_pangraph_differences(actual_pangraph, expected_pangraph)
            raise ex

    @staticmethod
    def show_pangraph_differences(actual_pangraph, expected_pangraph):
        print("Nodes differences: ")
        PangraphBuilderTests.show_nodes_differences(actual_pangraph.nodes,
                                                  expected_pangraph.nodes)

        print("Paths differences: ")
        PangraphBuilderTests.show_paths_differences(actual_pangraph.paths,
                                                  expected_pangraph.paths)

    @staticmethod
    def show_nodes_differences(actual_nodes, expected_nodes):
        if len(actual_nodes) != len(expected_nodes):
            print(f"Actual graph has {len(actual_nodes)} nodes while expected graph: {len(expected_nodes)}")
        for node_id in range(len(actual_nodes)):
            if actual_nodes[node_id] is None or actual_nodes[node_id] != expected_nodes[node_id]:
                print(f"Nodes {node_id} differ:")
                print(f"Actual: {actual_nodes[node_id]}")
                print(f"Expected: {expected_nodes[node_id]}")

    @staticmethod
    def show_paths_differences(actual_paths, expected_paths):
        for seq_id, paths in expected_paths.items():
            if seq_id not in actual_paths:
                print(f"There is no {seq_id} in actual paths!")
                continue
            elif actual_paths[seq_id] != expected_paths[seq_id]:
                print(f"Sequence {seq_id} differs in expected and actual paths!")
                print(f"Expected: {paths}")
                print(f"Actual: {actual_paths[seq_id]}")
        # if actual_paths.shape != expected_paths.shape:
        #     print(f"Actual pm has shape {actual_paths.shape} while expected pm: {expected_paths.shape}")
        #
        # for i, (actual_row, expected_row) in enumerate(zip(actual_paths, expected_paths)):
        #     if (actual_row != expected_row).any():
        #         print(f"Rows {i} differ.")
        #         print(f"Actual row: {actual_row}")
        #         print(f"Expected row: {expected_row}")
