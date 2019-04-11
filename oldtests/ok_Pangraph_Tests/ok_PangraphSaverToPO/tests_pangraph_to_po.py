import unittest
from ddt import ddt

from tests.Pangraph_Tests.Pangraph_Tests import PangraphTests

from tests.context import make_base as n
from tests.context import Node
from tests.context import MultialignmentMetadata
from tests.context import PangraphToPO
from tools import pathtools


@ddt
class PangraphSaverToPOTest_BuildPangraph(PangraphTests):

    def setUp(self):
        self.seq_metadata = MultialignmentMetadata(
            pathtools.get_file_content("Pangraph_Tests/seq_metadata.csv"))

    def test_1_typical_pangraph(self):
        expected_po_content_path = "Pangraph_Tests/" \
                                   "PangraphBuilderFromPO_Tests/" \
                                   "files_build_pangraph/" \
                                   "test_1.po"

        pangraph_nodes = [Node(node_id=0, base=n('A'), aligned_to=1),
                          Node(node_id=1, base=n('G'), aligned_to=0),
                          Node(node_id=2, base=n('C'), aligned_to=3),
                          Node(node_id=3, base=n('G'), aligned_to=2),
                          Node(node_id=4, base=n('A'), aligned_to=5),
                          Node(node_id=5, base=n('T'), aligned_to=4),
                          Node(node_id=6, base=n('G'), aligned_to=None),
                          Node(node_id=7, base=n('G'), aligned_to=None),
                          Node(node_id=8, base=n('A'), aligned_to=9),
                          Node(node_id=9, base=n('C'), aligned_to=10),
                          Node(node_id=10, base=n('G'), aligned_to=11),
                          Node(node_id=11, base=n('T'), aligned_to=8),
                          Node(node_id=12, base=n('A'), aligned_to=13),
                          Node(node_id=13, base=n('C'), aligned_to=12),
                          Node(node_id=14, base=n('T'), aligned_to=None),
                          Node(node_id=15, base=n('A'), aligned_to=16),
                          Node(node_id=16, base=n('C'), aligned_to=17),
                          Node(node_id=17, base=n('G'), aligned_to=15)
                         ]

        pangraph_paths = {
            'seq0': [[0, 2, 4, 6, 7, 8, 12, 14, 16]],
            'seq1': [[1, 2, 5, 6, 7, 9]],
            'seq2': [[3, 4, 6, 7, 10, 12, 14, 17]],
            'seq3': [[11, 13, 14, 15]]
        }

        pangraph = PangraphTests.setup_pangraph(pangraph_nodes, pangraph_paths)
        pangraph_to_po = PangraphToPO()
        actual_po_content = pangraph_to_po.get_po_file_content_from_pangraph(pangraph)

        expected_po_content = PangraphTests.get_file_content_as_str(expected_po_content_path)
        self.assertEqual(expected_po_content, actual_po_content)

    def test_2_consensuses_and_empty_sequences(self):
        expected_po_content_path = "Pangraph_Tests/" \
                                   "PangraphBuilderFromPO_Tests/" \
                                   "files_build_pangraph/" \
                                   "test_2.po"

        pangraph_nodes = [Node(node_id=0, base=n('C'), aligned_to=1),
                          Node(node_id=1, base=n('T'), aligned_to=0),
                          Node(node_id=2, base=n('A'), aligned_to=3),
                          Node(node_id=3, base=n('G'), aligned_to=2),
                          Node(node_id=4, base=n('C'), aligned_to=None),
                          Node(node_id=5, base=n('T'), aligned_to=None),
                          Node(node_id=6, base=n('A'), aligned_to=7),
                          Node(node_id=7, base=n('T'), aligned_to=6),
                          Node(node_id=8, base=n('G'), aligned_to=None)
                          ]

        pangraph_paths = {
            'seq0': [[0, 3, 4, 5, 6, 8]],
            'seq1': [[1, 2, 4, 5, 7, 8]],
            'seq2': [[]],
            'seq3': [[]],
            'CONSENS0': [[0, 3, 4, 5, 7, 8]],
            'CONSENS1': [[1, 2, 4, 5, 6, 8]]
        }
        pangraph = PangraphTests.setup_pangraph(pangraph_nodes, pangraph_paths)
        pangraph_to_po = PangraphToPO()
        actual_po_content = pangraph_to_po.get_po_file_content_from_pangraph(pangraph)

        expected_po_content = PangraphTests.get_file_content_as_str(expected_po_content_path)
        self.assertEqual(expected_po_content, actual_po_content)



if __name__ == '__main__':
    unittest.main()
