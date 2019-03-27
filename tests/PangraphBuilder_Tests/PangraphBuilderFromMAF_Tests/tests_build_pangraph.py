import unittest
from ddt import ddt

from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests

from tests.context import make_base as n
from tests.context import Node
from tests.context import MultialignmentMetadata
from tools import pathtools


@ddt
class PangraphBuilderFromMAFTest_BuildPangraph(PangraphBuilderTests):

    def setUp(self):
        self.seq_metadata = MultialignmentMetadata(
            pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.csv"))

    def test_1_messy_sequences(self):
        maf_path = "PangraphBuilder_Tests/" \
                   "PangraphBuilderFromMAF_Tests/" \
                   "files_build_pangraph/" \
                   "test_1_messy_sequences.maf"
        expected_nodes = [
            Node(node_id=0, base=n('A'), aligned_to=None),
            Node(node_id=1, base=n('A'), aligned_to=2),
            Node(node_id=2, base=n('C'), aligned_to=1),
            Node(node_id=3, base=n('T'), aligned_to=None),
            Node(node_id=4, base=n('C'), aligned_to=5),
            Node(node_id=5, base=n('G'), aligned_to=4),
            Node(node_id=6, base=n('A'), aligned_to=None),
            Node(node_id=7, base=n('C'), aligned_to=None),
            Node(node_id=8, base=n('G'), aligned_to=None),
            Node(node_id=9, base=n('C'), aligned_to=10),
            Node(node_id=10, base=n('G'), aligned_to=9),
            Node(node_id=11, base=n('T'), aligned_to=None),
            Node(node_id=12, base=n('C'), aligned_to=None),
            Node(node_id=13, base=n('C'), aligned_to=None),
            Node(node_id=14, base=n('A'), aligned_to=None),
        ]

        expected_paths = {
            "seq0": [[1, 3, 4, 6, 8, 9, 11, 12]],
            "seq1": [[2, 3, 4, 10, 11, 12, 13, 14]],
            "seq2": [[0, 2, 5, 6, 7, 10, 11, 12, 14]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf(maf_path, self.seq_metadata)
        self.assertEqual(expected_pangraph, actual_pangraph)


if __name__ == '__main__':
    unittest.main()
