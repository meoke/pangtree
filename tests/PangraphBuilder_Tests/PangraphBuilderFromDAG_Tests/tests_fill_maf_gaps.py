import unittest
from ddt import ddt

from tests.context import metadatareader
from tests.context import Node
from tests.context import nucleotides as n
from tests.context import pathtools
from tests.context import FastaSource

from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests


@ddt
class PangraphBuilderFromDAGTest_FillMafGaps(PangraphBuilderTests):
    class FakeFastaSource(FastaSource):
        def __init__(self):
            self.sources = {
                "seq0": "",
                "seq1": "ACTAGGT",
                "seq2": "GGTCAGT",
                "seq3": "",
                "seq4": ""}

        def get_source(self, id: str, start: int = None, end: int = None) -> str:
            try:
                return self.sources[id][start:end]
            except KeyError:
                raise Exception("No record found with given id!")

    @classmethod
    def setUpClass(cls):
        cls.seq_metadata = metadatareader.read(pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.json"))
        cls.fasta_source = PangraphBuilderFromDAGTest_FillMafGaps.FakeFastaSource()

    def test_1_missing_sequence_start(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_1_missing_sequence_start.maf"
        expected_nodes = [
            Node(id=0, base=n.code('A'), aligned_to=None),
            Node(id=1, base=n.code('C'), aligned_to=None),
            Node(id=2, base=n.code('T'), aligned_to=None),
            Node(id=3, base=n.code('A'), aligned_to=4),
            Node(id=4, base=n.code('G'), aligned_to=3),
            Node(id=5, base=n.code('G'), aligned_to=None),
            Node(id=6, base=n.code('G'), aligned_to=7),
            Node(id=7, base=n.code('T'), aligned_to=6),
            Node(id=8, base=n.code('C'), aligned_to=None),
            Node(id=9, base=n.code('A'), aligned_to=None),
            Node(id=10, base=n.code('G'), aligned_to=None),
            Node(id=11, base=n.code('T'), aligned_to=None)
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
            Node(id=0, base=n.code('A'), aligned_to=1),
            Node(id=1, base=n.code('G'), aligned_to=0),
            Node(id=2, base=n.code('C'), aligned_to=3),
            Node(id=3, base=n.code('G'), aligned_to=2),
            Node(id=4, base=n.code('T'), aligned_to=None),
            Node(id=5, base=n.code('A'), aligned_to=6),
            Node(id=6, base=n.code('C'), aligned_to=5),

            Node(id=7, base=n.code('A'), aligned_to=None),
            Node(id=8, base=n.code('G'), aligned_to=None),
            Node(id=9, base=n.code('G'), aligned_to=None),
            Node(id=10, base=n.code('T'), aligned_to=None),

            Node(id=11, base=n.code('G'), aligned_to=None),
            Node(id=12, base=n.code('T'), aligned_to=None),
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
            Node(id=0, base=n.code('A'), aligned_to=1),
            Node(id=1, base=n.code('G'), aligned_to=0),
            Node(id=2, base=n.code('C'), aligned_to=None),

            # missing seq1
            Node(id=3, base=n.code('T'), aligned_to=None),
            Node(id=4, base=n.code('A'), aligned_to=None),

            # missing seq2
            Node(id=5, base=n.code('G'), aligned_to=None),
            Node(id=6, base=n.code('T'), aligned_to=None),

            # block 1
            Node(id=7, base=n.code('C'), aligned_to=8),
            Node(id=8, base=n.code('G'), aligned_to=7),
            Node(id=9, base=n.code('A'), aligned_to=None),
            Node(id=10, base=n.code('G'), aligned_to=None),
            Node(id=11, base=n.code('T'), aligned_to=None)
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
            Node(id=0, base=n.code('A'), aligned_to=1),
            Node(id=1, base=n.code('G'), aligned_to=0),
            Node(id=2, base=n.code('C'), aligned_to=None),
            Node(id=3, base=n.code('T'), aligned_to=None),
            Node(id=4, base=n.code('A'), aligned_to=5),
            Node(id=5, base=n.code('G'), aligned_to=4),

            # missing se2
            Node(id=6, base=n.code('T'), aligned_to=None),
            Node(id=7, base=n.code('C'), aligned_to=None),

            Node(id=8, base=n.code('A'), aligned_to=9),
            Node(id=9, base=n.code('G'), aligned_to=8),
            Node(id=10, base=n.code('G'), aligned_to=None),
            Node(id=11, base=n.code('T'), aligned_to=None),


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

    def test_5_missing_one_reverted_sequence_middle(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/" \
                   "test_5_missing_one_reverted_sequence_middle.maf"
        expected_nodes = [
            # block 0
            Node(id=0, base=n.code('A'), aligned_to=1),
            Node(id=1, base=n.code('G'), aligned_to=0),
            Node(id=2, base=n.code('C'), aligned_to=None),
            Node(id=3, base=n.code('T'), aligned_to=None),
            Node(id=4, base=n.code('A'), aligned_to=5),
            Node(id=5, base=n.code('G'), aligned_to=4),

            # missing seq2
            # Node(id=6, base=n.code('T'), in_nodes=[5], aligned_to=None),
            # Node(id=7, base=n.code('C'), in_nodes=[10], aligned_to=None),

            Node(id=6, base=n.code('A'), aligned_to=7),
            Node(id=7, base=n.code('G'), aligned_to=6),
            Node(id=8, base=n.code('G'), aligned_to=None),
            Node(id=9, base=n.code('T'), aligned_to=None),
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 2, 3, 4, 7, 8, 9]],
            "seq2": [[1, 5], [6, 8, 9]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path,
                                                                                                self.seq_metadata,
                                                                                                self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


if __name__ == '__main__':
    unittest.main()
