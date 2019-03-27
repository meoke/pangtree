import unittest
from ddt import ddt

from tests.context import MultialignmentMetadata
from tests.context import Node
from tests.context import make_nucleobase as n
from tests.context import pathtools
from tests.context import FastaProvider

from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests


@ddt
class PangraphBuilderFromDAGTest_FillMafGapsWithN(PangraphBuilderTests):
    class FakeFastaSource(FastaProvider):
        def __init__(self):
            self.sources = {
                "seq0": "",
                "seq1": "ACTAGGT",
                "seq2": "GGTCAGT",
                "seq3": "",
                "seq4": ""}

        def get_source(self, sequence_id: str, start: int = None, end: int = None) -> str:
            try:
                return self.sources[sequence_id][start:end]
            except KeyError:
                raise Exception("No record found with given id!")

    @classmethod
    def setUpClass(cls):
        cls.unknown_symbol = '?'
        cls.seq_metadata = MultialignmentMetadata(pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.csv"))
        cls.fasta_source = None

    def test_1_missing_sequence_start(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_1_missing_sequence_start.maf"
        expected_nodes = [
            Node(node_id=0, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=1, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=2, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=3, base=n('A'), aligned_to=4),
            Node(node_id=4, base=n('G'), aligned_to=3),
            Node(node_id=5, base=n('G'), aligned_to=None),
            Node(node_id=6, base=n('G'), aligned_to=7),
            Node(node_id=7, base=n('T'), aligned_to=6),
            Node(node_id=8, base=n('C'), aligned_to=None),
            Node(node_id=9, base=n('A'), aligned_to=None),
            Node(node_id=10, base=n('G'), aligned_to=None),
            Node(node_id=11, base=n('T'), aligned_to=None)
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 1, 2, 3, 5, 6, 11]],
            "seq2": [[4, 5, 7, 8, 9, 10, 11]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_2_missing_sequence_end(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_2_missing_sequence_end.maf"
        expected_nodes = [
            Node(node_id=0, base=n('A'), aligned_to=1),
            Node(node_id=1, base=n('G'), aligned_to=0),
            Node(node_id=2, base=n('C'), aligned_to=3),
            Node(node_id=3, base=n('G'), aligned_to=2),
            Node(node_id=4, base=n('T'), aligned_to=None),
            Node(node_id=5, base=n('A'), aligned_to=6),
            Node(node_id=6, base=n('C'), aligned_to=5),

            Node(node_id=7, base=n('A'), aligned_to=None),
            Node(node_id=8, base=n('G'), aligned_to=None),
            Node(node_id=9, base=n('G'), aligned_to=None),
            Node(node_id=10, base=n('T'), aligned_to=None),

            Node(node_id=11, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=12, base=n(self.unknown_symbol), aligned_to=None),
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 2, 4, 5, 8, 9, 10]],
            "seq2": [[1, 3, 4, 6, 7, 11, 12]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_3_missing_two_sequences_middle(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_3_missing_two_sequences_middle.maf"
        expected_nodes = [
            # block 0
            Node(node_id=0, base=n('A'), aligned_to=1),
            Node(node_id=1, base=n('G'), aligned_to=0),
            Node(node_id=2, base=n('C'), aligned_to=None),

            # missing seq1
            Node(node_id=3, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=4, base=n(self.unknown_symbol), aligned_to=None),

            # missing seq2
            Node(node_id=5, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=6, base=n(self.unknown_symbol), aligned_to=None),

            # block 1
            Node(node_id=7, base=n('C'), aligned_to=8),
            Node(node_id=8, base=n('G'), aligned_to=7),
            Node(node_id=9, base=n('A'), aligned_to=None),
            Node(node_id=10, base=n('G'), aligned_to=None),
            Node(node_id=11, base=n('T'), aligned_to=None)
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 2, 3, 4, 8, 10, 11]],
            "seq2": [[1, 5, 6, 7, 9, 10, 11]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_4_missing_one_sequence_middle(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/" \
                   "test_4_missing_one_sequence_middle.maf"
        expected_nodes = [
            #block 0
            Node(node_id=0, base=n('A'), aligned_to=1),
            Node(node_id=1, base=n('G'), aligned_to=0),
            Node(node_id=2, base=n('C'), aligned_to=None),
            Node(node_id=3, base=n('T'), aligned_to=None),
            Node(node_id=4, base=n('A'), aligned_to=5),
            Node(node_id=5, base=n('G'), aligned_to=4),

            # missing se2
            Node(node_id=6, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=7, base=n(self.unknown_symbol), aligned_to=None),

            Node(node_id=8, base=n('A'), aligned_to=9),
            Node(node_id=9, base=n('G'), aligned_to=8),
            Node(node_id=10, base=n('G'), aligned_to=None),
            Node(node_id=11, base=n('T'), aligned_to=None),


        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 2, 3, 4, 9, 10, 11]],
            "seq2": [[1, 5, 6, 7, 8, 10, 11]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_5_missing_one_reverted_sequence_middle_1_1(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/" \
                   "test_5_missing_one_reverted_sequence_middle_1_1.maf"
        expected_nodes = [
            # block 0
            Node(node_id=0, base=n('A'), aligned_to=1),
            Node(node_id=1, base=n('G'), aligned_to=0),
            Node(node_id=2, base=n('C'), aligned_to=None),
            Node(node_id=3, base=n('T'), aligned_to=None),
            Node(node_id=4, base=n('A'), aligned_to=5),
            Node(node_id=5, base=n('G'), aligned_to=4),

            # missing seq2, on edge (1,1)
            Node(node_id=6, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=7, base=n(self.unknown_symbol), aligned_to=None),

            # block 1
            Node(node_id=8, base=n('A'), aligned_to=9),
            Node(node_id=9, base=n('G'), aligned_to=8),
            Node(node_id=10, base=n('G'), aligned_to=None),
            Node(node_id=11, base=n('T'), aligned_to=None),
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 2, 3, 4, 9, 10, 11]],
            "seq2": [[1, 5, 6, 7], [8, 10, 11]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path,
                                                                                                self.seq_metadata,
                                                                                                self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_6_missing_one_reverted_sequence_middle_minus1_1(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/" \
                   "test_6_missing_one_reverted_sequence_middle_minus1_1.maf"

        expected_nodes = [
            # block 1 because it is first in DAG and reverted
            Node(node_id=0, base=n('A'), aligned_to=None),
            Node(node_id=1, base=n('C'), aligned_to=None),
            Node(node_id=2, base=n('C'), aligned_to=3),
            Node(node_id=3, base=n('T'), aligned_to=2),

            # missing seq2, on edge (-1,1)
            Node(node_id=4, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=5, base=n(self.unknown_symbol), aligned_to=None),

            Node(node_id=6, base=n('A'), aligned_to=7),
            Node(node_id=7, base=n('C'), aligned_to=6),
            Node(node_id=8, base=n('C'), aligned_to=None),
            Node(node_id=9, base=n('T'), aligned_to=None),
            Node(node_id=10, base=n('A'), aligned_to=11),
            Node(node_id=11, base=n('C'), aligned_to=10),


        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 1, 2], [6, 8, 9, 10]],
            "seq2": [[0, 1, 3, 4, 5, 7, 11]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path,
                                                                                                self.seq_metadata,
                                                                                                self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_7_missing_one_reverted_sequence_middle_minus1_minus1(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/" \
                   "test_7_missing_one_reverted_sequence_middle_minus1_minus1.maf"
        expected_nodes = [
            # block 0
            Node(node_id=0, base=n('A'), aligned_to=None),
            Node(node_id=1, base=n('C'), aligned_to=None),
            Node(node_id=2, base=n('T'), aligned_to=None),
            Node(node_id=3, base=n('A'), aligned_to=None),

            # missing seq2
            Node(node_id=4, base=n(self.unknown_symbol), aligned_to=None),
            Node(node_id=5, base=n(self.unknown_symbol), aligned_to=None),

            # block 1
            Node(node_id=6, base=n('A'), aligned_to=7),
            Node(node_id=7, base=n('G'), aligned_to=6),
            Node(node_id=8, base=n('C'), aligned_to=9),
            Node(node_id=9, base=n('G'), aligned_to=8),
            Node(node_id=10, base=n('C'), aligned_to=11),
            Node(node_id=11, base=n('T'), aligned_to=10),
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 1, 2, 3, 7, 9, 11]],
            "seq2": [[0, 1, 4, 5, 6, 8, 10]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path,
                                                                                                self.seq_metadata,
                                                                                                self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


if __name__ == '__main__':
    unittest.main()
