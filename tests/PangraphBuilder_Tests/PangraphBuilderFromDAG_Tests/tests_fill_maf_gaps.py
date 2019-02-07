import unittest
from ddt import ddt

from context import metadatareader
from context import Node
from context import nucleotides as n
from context import pathtools
from context import FastaSource

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
                self.sources[id][start:end]
            except KeyError:
                raise Exception("No record found with given id!")

    @classmethod
    def setUpClass(cls):
        cls.seq_metadata = metadatareader.read(pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.json"))
        cls.fasta_source = PangraphBuilderFromDAGTest_FillMafGaps.FakeFastaSource()

    def test_1_missing_sequence_start(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_1_missing_sequence_start.maf"
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

        expected_paths = {
            "testseq1": [0, 1, 2, 3, 5, 6, 11],
            "testseq2": [4, 5, 7, 8, 9, 10, 11]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_2_missing_sequence_end(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_2_missing_sequence_end.maf"
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=1),
            Node(id=1, base=n.code('G'), in_nodes=[], aligned_to=0),
            Node(id=2, base=n.code('C'), in_nodes=[0], aligned_to=3),
            Node(id=3, base=n.code('G'), in_nodes=[1], aligned_to=2),
            Node(id=4, base=n.code('T'), in_nodes=[2,3], aligned_to=None),
            Node(id=5, base=n.code('A'), in_nodes=[4], aligned_to=6),
            Node(id=6, base=n.code('C'), in_nodes=[4], aligned_to=5),
            Node(id=7, base=n.code('A'), in_nodes=[6], aligned_to=None),
            Node(id=8, base=n.code('G'), in_nodes=[5], aligned_to=None),
            Node(id=9, base=n.code('G'), in_nodes=[8], aligned_to=None),
            Node(id=10, base=n.code('T'), in_nodes=[9], aligned_to=None),
            Node(id=11, base=n.code('G'), in_nodes=[7], aligned_to=None),
            Node(id=12, base=n.code('T'), in_nodes=[11], aligned_to=None),
        ]

        expected_paths = {
            "testseq1": [0, 2, 4, 5, 8, 9, 10],
            "testseq2": [1, 3, 4, 6, 7, 11, 12]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_3_missing_two_sequences_middle(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_3_missing_two_sequences_middle.maf"
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=1),
            Node(id=1, base=n.code('G'), in_nodes=[], aligned_to=0),
            Node(id=2, base=n.code('C'), in_nodes=[0], aligned_to=None),

            Node(id=3, base=n.code('C'), in_nodes=[11], aligned_to=4),
            Node(id=4, base=n.code('G'), in_nodes=[9], aligned_to=3),
            Node(id=5, base=n.code('A'), in_nodes=[3], aligned_to=None),
            Node(id=6, base=n.code('G'), in_nodes=[4, 5], aligned_to=None),
            Node(id=7, base=n.code('T'), in_nodes=[6], aligned_to=None),

            Node(id=8, base=n.code('T'), in_nodes=[2], aligned_to=None),
            Node(id=9, base=n.code('A'), in_nodes=[8], aligned_to=None),

            Node(id=10, base=n.code('G'), in_nodes=[1], aligned_to=None),
            Node(id=11, base=n.code('T'), in_nodes=[10], aligned_to=None)
        ]

        expected_paths = {
            "testseq1": [*sorted([0, 2, 8, 9, 4, 6, 7])],
            "testseq2": [*sorted([1, 10, 11, 3, 5, 6, 7])]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


    def test_4_missing_one_sequence_middle(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/files_fill_maf_gaps/test_4_missing_one_sequence_middle"
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=1),
            Node(id=1, base=n.code('G'), in_nodes=[], aligned_to=0),
            Node(id=2, base=n.code('C'), in_nodes=[0], aligned_to=None),
            Node(id=3, base=n.code('T'), in_nodes=[2], aligned_to=None),
            Node(id=4, base=n.code('A'), in_nodes=[3], aligned_to=5),
            Node(id=5, base=n.code('G'), in_nodes=[1], aligned_to=4),

            Node(id=6, base=n.code('A'), in_nodes=[11], aligned_to=7),
            Node(id=7, base=n.code('G'), in_nodes=[4], aligned_to=6),
            Node(id=8, base=n.code('G'), in_nodes=[6, 7], aligned_to=None),
            Node(id=9, base=n.code('T'), in_nodes=[8], aligned_to=None),

            Node(id=10, base=n.code('T'), in_nodes=[5], aligned_to=None),
            Node(id=11, base=n.code('C'), in_nodes=[10], aligned_to=None)
        ]

        expected_paths = {
            "testseq1": [*sorted([0, 2, 3, 4, 7, 8, 9])],
            "testseq2": [*sorted([1, 5, 10, 11, 6, 8, 9])]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


if __name__ == '__main__':
    unittest.main()
