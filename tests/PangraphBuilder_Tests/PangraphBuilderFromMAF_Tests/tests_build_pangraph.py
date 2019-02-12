import unittest
from ddt import ddt

from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests

from tests.context import nucleotides as n
from tests.context import Node
from tests.context import metadatareader
from tests.context import pathtools

@ddt
class PangraphBuilderFromMAFTest_BuildPangraph(PangraphBuilderTests):

    def setUp(self):
        self.seq_metadata = metadatareader.read(
            pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.json"))

    def test_1_messy_sequences(self):
        maf_path = "PangraphBuilder_Tests/" \
                   "PangraphBuilderFromMAF_Tests/" \
                   "files_build_pangraph/" \
                   "test_1_messy_sequences.maf"
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('A'), in_nodes=[], aligned_to=2),
            Node(id=2, base=n.code('C'), in_nodes=[0], aligned_to=1),
            Node(id=3, base=n.code('T'), in_nodes=[1, 2], aligned_to=None),
            Node(id=4, base=n.code('C'), in_nodes=[3], aligned_to=5),
            Node(id=5, base=n.code('G'), in_nodes=[2], aligned_to=4),
            Node(id=6, base=n.code('A'), in_nodes=[4, 5], aligned_to=None),
            Node(id=7, base=n.code('C'), in_nodes=[6], aligned_to=None),
            Node(id=8, base=n.code('G'), in_nodes=[6], aligned_to=None),
            Node(id=9, base=n.code('C'), in_nodes=[8], aligned_to=10),
            Node(id=10, base=n.code('G'), in_nodes=[4, 7], aligned_to=9),
            Node(id=11, base=n.code('T'), in_nodes=[9, 10], aligned_to=None),
            Node(id=12, base=n.code('C'), in_nodes=[11], aligned_to=None),
            Node(id=13, base=n.code('C'), in_nodes=[12], aligned_to=None),
            Node(id=14, base=n.code('A'), in_nodes=[12, 13], aligned_to=None),
        ]

        expected_paths = {
            "seq0": [1, 3, 4, 6, 8, 9, 11, 12],
            "seq1": [2, 3, 4, 10, 11, 12, 13, 14],
            "seq2": [0, 2, 5, 6, 7, 10, 11, 12, 14]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf(maf_path, self.seq_metadata)
        self.assertEqual(expected_pangraph, actual_pangraph)


if __name__ == '__main__':
    unittest.main()
