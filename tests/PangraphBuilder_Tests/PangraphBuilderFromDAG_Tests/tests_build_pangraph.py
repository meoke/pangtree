import unittest
from ddt import ddt

from context import metadatareader
from context import Node
from context import nucleotides as n
from context import pathtools

from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests


@ddt
class PangraphBuilderFromDAGTest_BuildPangraph(PangraphBuilderTests):
    @classmethod
    def setUpClass(cls):
        cls.seq_metadata = metadatareader.read(
            pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.json"))
        cls.fasta_source = None

    def test_0_simple(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_0_simple.maf"
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

        expected_paths = {
            "seq0": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            "seq1": [0, 1, 2, 3, 4, 5, 6, 7, 8]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_1_reversed_seq_in_one_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_1_reversed_seq_in_one_block.maf"
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=1),
            Node(id=1, base=n.code('G'), in_nodes=[], aligned_to=0),
            Node(id=2, base=n.code('C'), in_nodes=[0], aligned_to=3),
            Node(id=3, base=n.code('T'), in_nodes=[1], aligned_to=2),

            Node(id=4, base=n.code('G'), in_nodes=[2], aligned_to=5),
            Node(id=5, base=n.code('T'), in_nodes=[], aligned_to=4),
            Node(id=6, base=n.code('A'), in_nodes=[5], aligned_to=7),
            Node(id=7, base=n.code('C'), in_nodes=[4], aligned_to=6),
            Node(id=8, base=n.code('C'), in_nodes=[7], aligned_to=None),

        ]

        expected_pats = {
            "seq0": [0, 2, 4, 7, 8],
            "seq1": [1, 3, 5, 6]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_pats)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_2_seq_starts_in_second_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_2_seq_starts_in_second_block.maf"
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

        expected_paths = {
            "seq0": [1, 2, 3],
            "seq1": [0, 5, 7],
            "seq2": [3, 4, 6, 8]
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_3_edge_not_from_last_node_in_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_3_edge_not_from_last_node_in_block.maf"
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
            "seq1": [0, 1, 3, 4, 6, 7],
            "seq2": [0, 1, 2, 3, 5],
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_pats)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_4_single_block_no_nucleotides(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_4_single_block_no_nucleotides.maf"
        expected_nodes = []

        expected_pats = {
            "seq0": [],
            "seq1": []
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_pats)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_5_single_block_single_nucletodide(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_5_single_block_single_nucletodide.maf"
        expected_nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None)]

        expected_paths = {
            "seq0": [0],
            "seq1": [0],
            "seq2": [0],
            "seq3": [0]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_6_1st_block_separates_into_2_branches_which_connect_in_3rd_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/" \
                   "test_6_1st_block_separates_into_2_branches_which_connect_in_3rd_block.maf"
        expected_nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=1),
                          Node(id=1, base=n.code('C'), in_nodes=[], aligned_to=2),
                          Node(id=2, base=n.code('G'), in_nodes=[], aligned_to=0),
                          Node(id=3, base=n.code('C'), in_nodes=[0, 1, 2], aligned_to=None),
                          Node(id=4, base=n.code('A'), in_nodes=[3], aligned_to=5),
                          Node(id=5, base=n.code('T'), in_nodes=[3], aligned_to=4),
                          Node(id=6, base=n.code('C'), in_nodes=[4], aligned_to=7),
                          Node(id=7, base=n.code('G'), in_nodes=[10], aligned_to=8),
                          Node(id=8, base=n.code('T'), in_nodes=[5], aligned_to=6),
                          Node(id=9, base=n.code('G'), in_nodes=[5], aligned_to=None),
                          Node(id=10, base=n.code('G'), in_nodes=[9], aligned_to=None)]

        expected_paths = {
            "seq0": [0, 3, 4, 6],
            "seq1": [1, 3, 5, 7, 9, 10],
            "seq2": [2, 3, 5, 8]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_7_inactive_edges_due_to_reversed_seqs(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_7_inactive_edges_due_to_reversed_seqs.maf"
        expected_nodes = [
                            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
                            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
                            Node(id=2, base=n.code('T'), in_nodes=[0, 1], aligned_to=None),
                            Node(id=3, base=n.code('G'), in_nodes=[], aligned_to=None),
                            Node(id=4, base=n.code('A'), in_nodes=[2], aligned_to=5),
                            Node(id=5, base=n.code('G'), in_nodes=[3], aligned_to=4),
                            Node(id=6, base=n.code('C'), in_nodes=[4], aligned_to=None),
                            Node(id=7, base=n.code('A'), in_nodes=[6], aligned_to=None),
                            Node(id=8, base=n.code('T'), in_nodes=[6, 7], aligned_to=None),
                            Node(id=9, base=n.code('G'), in_nodes=[8], aligned_to=None),
                            Node(id=10, base=n.code('A'), in_nodes=[9], aligned_to=None),
                            Node(id=11, base=n.code('A'), in_nodes=[9, 10], aligned_to=None),
                            Node(id=12, base=n.code('A'), in_nodes=[11], aligned_to=13),
                            Node(id=13, base=n.code('C'), in_nodes=[11], aligned_to=12),
                            Node(id=14, base=n.code('A'), in_nodes=[12], aligned_to=15),
                            Node(id=15, base=n.code('T'), in_nodes=[12, 13], aligned_to=14),
        ]

        expected_paths = {
            "seq1": [0, 1, 2, 4, 6, 8, 9, 10, 11, 12, 14],
            "seq2": [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 13, 15],
            "seq3": [0, 2, 3, 5, 6, 7, 8, 9, 11, 12, 15]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_8_reversed_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_8_reversed_block.maf"
        expected_nodes = [
                            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
                            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
                            Node(id=2, base=n.code('T'), in_nodes=[0,1], aligned_to=None),
                            #next block is reversed because it was converted to dag
                            Node(id=3, base=n.code('C'), in_nodes=[2], aligned_to=4),
                            Node(id=4, base=n.code('T'), in_nodes=[2], aligned_to=3),
                            Node(id=5, base=n.code('C'), in_nodes=[3,4], aligned_to=None),
                            Node(id=6, base=n.code('C'), in_nodes=[5], aligned_to=None),
                            Node(id=7, base=n.code('A'), in_nodes=[6], aligned_to=None),
                            Node(id=8, base=n.code('T'), in_nodes=[6,7], aligned_to=None),
                            Node(id=9, base=n.code('G'), in_nodes=[8], aligned_to=None)
        ]

        expected_paths = {
            "seq1": [0, 1, 2, 4, 5, 6, 8, 9],
            "seq2": [0, 1, 2, 3, 5, 6, 7, 8, 9],
            "seq3": [0, 2, 3, 5, 6, 7, 8, 9]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_9_inactive_edges_but_all_strands_plus(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/" \
                   "test_9_inactive_edges_but_all_strands_plus.maf"
        expected_nodes = [
            Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
            Node(id=1, base=n.code('C'), in_nodes=[0], aligned_to=None),
            Node(id=2, base=n.code('T'), in_nodes=[1], aligned_to=None),
            Node(id=3, base=n.code('G'), in_nodes=[2], aligned_to=None),
            Node(id=4, base=n.code('G'), in_nodes=[3], aligned_to=None),

            Node(id=5, base=n.code('A'), in_nodes=[4], aligned_to=None),
            Node(id=6, base=n.code('C'), in_nodes=[5], aligned_to=None),
            Node(id=7, base=n.code('T'), in_nodes=[6], aligned_to=None),
            Node(id=8, base=n.code('G'), in_nodes=[7], aligned_to=None),
            Node(id=9, base=n.code('G'), in_nodes=[8], aligned_to=None),

            Node(id=10, base=n.code('A'), in_nodes=[4, 9], aligned_to=None),
            Node(id=11, base=n.code('C'), in_nodes=[10], aligned_to=None),
            Node(id=12, base=n.code('T'), in_nodes=[11], aligned_to=None),
            Node(id=13, base=n.code('G'), in_nodes=[12], aligned_to=None),
            Node(id=14, base=n.code('G'), in_nodes=[13], aligned_to=None),

            Node(id=15, base=n.code('A'), in_nodes=[9, 14], aligned_to=None),
            Node(id=16, base=n.code('C'), in_nodes=[15], aligned_to=None),
            Node(id=17, base=n.code('T'), in_nodes=[16], aligned_to=None),
            Node(id=18, base=n.code('G'), in_nodes=[17], aligned_to=None),
            Node(id=19, base=n.code('G'), in_nodes=[18], aligned_to=None)
        ]

        expected_paths = {
            "seq1": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
            "seq2": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


if __name__ == '__main__':
    unittest.main()
