import unittest
from ddt import ddt

from tests.context import metadatareader
from tests.context import Node
from tests.context import nucleotides as n
from tests.context import pathtools

from tests.PangraphBuilder_Tests.PangraphBuilder_Tests import PangraphBuilderTests


@ddt
class PangraphBuilderFromDAGTest_BuildPangraph(PangraphBuilderTests):
    @classmethod
    def setUpClass(cls):
        cls.seq_metadata = metadatareader.read(
            pathtools.get_file_content("PangraphBuilder_Tests/seq_metadata.json"))
        cls.fasta_source = None

    def test_00_simple(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_0_simple.maf"
        expected_nodes = [
            Node(node_id=0, base=n.code('A'), aligned_to=None),
            Node(node_id=1, base=n.code('C'), aligned_to=None),
            Node(node_id=2, base=n.code('T'), aligned_to=None),
            Node(node_id=3, base=n.code('G'), aligned_to=None),
            Node(node_id=4, base=n.code('A'), aligned_to=None),
            Node(node_id=5, base=n.code('C'), aligned_to=None),
            Node(node_id=6, base=n.code('T'), aligned_to=None),
            Node(node_id=7, base=n.code('G'), aligned_to=None),
            Node(node_id=8, base=n.code('A'), aligned_to=None),
            Node(node_id=9, base=n.code('A'), aligned_to=None),
        ]

        expected_paths = {
            "seq0": [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]],
            "seq1": [[0, 1, 2, 3, 4, 5, 6, 7, 8]],
            "seq2": [],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_01_reversed_seq_in_one_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_1_reversed_seq_in_one_block.maf"
        expected_nodes = [
            Node(node_id=0, base=n.code('A'), aligned_to=1),
            Node(node_id=1, base=n.code('G'), aligned_to=0),
            Node(node_id=2, base=n.code('C'), aligned_to=3),
            Node(node_id=3, base=n.code('T'), aligned_to=2),

            Node(node_id=4, base=n.code('G'), aligned_to=5),
            Node(node_id=5, base=n.code('T'), aligned_to=4),
            Node(node_id=6, base=n.code('A'), aligned_to=7),
            Node(node_id=7, base=n.code('C'), aligned_to=6),
            Node(node_id=8, base=n.code('C'), aligned_to=None),

        ]
        expected_paths = {
            "seq0": [[0, 2, 4, 7, 8]],
            "seq1": [[1, 3], [5, 6]],
            "seq2": [],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_02_seq_starts_in_second_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_2_seq_starts_in_second_block.maf"
        expected_nodes = [
            Node(node_id=0, base=n.code('C'), aligned_to=None),
            Node(node_id=1, base=n.code('T'), aligned_to=None),
            Node(node_id=2, base=n.code('G'), aligned_to=None),
            Node(node_id=3, base=n.code('T'), aligned_to=None),
            Node(node_id=4, base=n.code('G'), aligned_to=5),
            Node(node_id=5, base=n.code('T'), aligned_to=4),
            Node(node_id=6, base=n.code('A'), aligned_to=None),
            Node(node_id=7, base=n.code('A'), aligned_to=None),
            Node(node_id=8, base=n.code('C'), aligned_to=None),
        ]

        expected_paths = {
            "seq0": [[1, 2, 3]],
            "seq1": [[0, 5, 7]],
            "seq2": [[3, 4, 6, 8]],
            "seq3": []
        }
        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_03_edge_not_from_last_node_in_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_3_edge_not_from_last_node_in_block.maf"
        expected_nodes = [
            Node(node_id=0, base=n.code('A'), aligned_to=None),
            Node(node_id=1, base=n.code('C'), aligned_to=None),
            Node(node_id=2, base=n.code('T'), aligned_to=None),
            Node(node_id=3, base=n.code('G'), aligned_to=None),
            Node(node_id=4, base=n.code('G'), aligned_to=None),
            Node(node_id=5, base=n.code('A'), aligned_to=6),
            Node(node_id=6, base=n.code('T'), aligned_to=5),
            Node(node_id=7, base=n.code('C'), aligned_to=None)
        ]

        expected_pats = {
            "seq1": [[0, 1, 3, 4, 6, 7]],
            "seq2": [[0, 1, 2, 3, 5]],
            "seq0": [],
            "seq3": []
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_pats)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_04_single_block_no_nucleotides(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_4_single_block_no_nucleotides.maf"
        expected_nodes = []

        expected_paths = {
            "seq0": [],
            "seq1": [],
            "seq2": [],
            "seq3": []
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_05_single_block_single_nucletodide(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_5_single_block_single_nucletodide.maf"
        expected_nodes = [Node(node_id=0, base=n.code('A'), aligned_to=None)]

        expected_paths = {
            "seq0": [[0]],
            "seq1": [[0]],
            "seq2": [[0]],
            "seq3": [[0]]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_06_1st_block_separates_into_2_branches_which_connect_in_3rd_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/" \
                   "test_6_1st_block_separates_into_2_branches_which_connect_in_3rd_block.maf"
        expected_nodes = [Node(node_id=0, base=n.code('A'), aligned_to=1),
                          Node(node_id=1, base=n.code('C'), aligned_to=2),
                          Node(node_id=2, base=n.code('G'), aligned_to=0),
                          Node(node_id=3, base=n.code('C'), aligned_to=None),
                          Node(node_id=4, base=n.code('A'), aligned_to=5),
                          Node(node_id=5, base=n.code('T'), aligned_to=4),

                          Node(node_id=6, base=n.code('G'), aligned_to=None),
                          Node(node_id=7, base=n.code('G'), aligned_to=None),

                          Node(node_id=8, base=n.code('C'), aligned_to=9),
                          Node(node_id=9, base=n.code('G'), aligned_to=10),
                          Node(node_id=10, base=n.code('T'), aligned_to=8)]

        expected_paths = {
            "seq0": [[0, 3, 4, 8]],
            "seq1": [[1, 3, 5, 6, 7, 9]],
            "seq2": [[2, 3, 5, 10]],
            "seq3": []
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_07_inactive_edges_due_to_reversed_seqs(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_7_inactive_edges_due_to_reversed_seqs.maf"
        expected_nodes = [
                            Node(node_id=0, base=n.code('A'), aligned_to=None),
                            Node(node_id=1, base=n.code('C'), aligned_to=None),
                            Node(node_id=2, base=n.code('T'), aligned_to=None),
                            Node(node_id=3, base=n.code('G'), aligned_to=None),
                            Node(node_id=4, base=n.code('A'), aligned_to=5),
                            Node(node_id=5, base=n.code('G'), aligned_to=4),
                            Node(node_id=6, base=n.code('C'), aligned_to=None),
                            Node(node_id=7, base=n.code('A'), aligned_to=None),
                            Node(node_id=8, base=n.code('T'), aligned_to=None),
                            Node(node_id=9, base=n.code('G'), aligned_to=None),
                            Node(node_id=10, base=n.code('A'), aligned_to=None),
                            Node(node_id=11, base=n.code('A'), aligned_to=None),
                            Node(node_id=12, base=n.code('A'), aligned_to=13),
                            Node(node_id=13, base=n.code('C'), aligned_to=12),
                            Node(node_id=14, base=n.code('A'), aligned_to=15),
                            Node(node_id=15, base=n.code('T'), aligned_to=14),
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 1, 2, 4, 6, 8, 9, 10, 11, 12, 14]],
            "seq2": [[0, 1, 2], [3, 5], [6, 7, 8, 9, 10, 11, 13, 15]],
            "seq3": [[0, 2], [3, 5], [6, 7, 8, 9, 11, 12, 15]]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_08_reversed_block(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/test_8_reversed_block.maf"
        expected_nodes = [
                            Node(node_id=0, base=n.code('C'), aligned_to=None),
                            Node(node_id=1, base=n.code('A'), aligned_to=None),
                            Node(node_id=2, base=n.code('T'), aligned_to=None),
                            #next block is reversed because it was converted to dag
                            Node(node_id=3, base=n.code('G'), aligned_to=None),
                            Node(node_id=4, base=n.code('G'), aligned_to=None),
                            Node(node_id=5, base=n.code('A'), aligned_to=6),
                            Node(node_id=6, base=n.code('G'), aligned_to=5),
                            Node(node_id=7, base=n.code('A'), aligned_to=None),
                            Node(node_id=8, base=n.code('G'), aligned_to=None),
                            Node(node_id=9, base=n.code('T'), aligned_to=None)
        ]

        expected_paths = {
            "seq0": [],
            "seq1": [[0, 1, 3, 4, 5, 7, 8, 9]],
            "seq2": [[0, 1, 2, 3, 4, 6, 7, 8, 9]],
            "seq3": [[0, 1, 2, 3, 4, 6, 7, 9]]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_09_inactive_edges_but_all_strands_plus(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/" \
                   "test_9_inactive_edges_but_all_strands_plus.maf"
        expected_nodes = [
            Node(node_id=0, base=n.code('A'), aligned_to=None),
            Node(node_id=1, base=n.code('C'), aligned_to=None),
            Node(node_id=2, base=n.code('T'), aligned_to=None),
            Node(node_id=3, base=n.code('G'), aligned_to=None),
            Node(node_id=4, base=n.code('G'), aligned_to=None),

            Node(node_id=5, base=n.code('A'), aligned_to=None),
            Node(node_id=6, base=n.code('C'), aligned_to=None),
            Node(node_id=7, base=n.code('T'), aligned_to=None),
            Node(node_id=8, base=n.code('G'), aligned_to=None),
            Node(node_id=9, base=n.code('G'), aligned_to=None),

            Node(node_id=10, base=n.code('A'), aligned_to=None),
            Node(node_id=11, base=n.code('C'), aligned_to=None),
            Node(node_id=12, base=n.code('T'), aligned_to=None),
            Node(node_id=13, base=n.code('G'), aligned_to=None),
            Node(node_id=14, base=n.code('G'), aligned_to=None),

            Node(node_id=15, base=n.code('A'), aligned_to=None),
            Node(node_id=16, base=n.code('C'), aligned_to=None),
            Node(node_id=17, base=n.code('T'), aligned_to=None),
            Node(node_id=18, base=n.code('G'), aligned_to=None),
            Node(node_id=19, base=n.code('G'), aligned_to=None)
        ]
        expected_paths = {
            "seq0": [],
            "seq1": [[0, 1, 2, 3, 4, 10, 11, 12, 13, 14], [5, 6, 7, 8, 9, 15, 16, 17, 18, 19]],
            "seq2": [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]],
            "seq3": []
        }


        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path, self.seq_metadata, self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)

    def test_10_parallel_blocks_1st_and_2nd_merge_into_3rd(self):
        maf_path = "PangraphBuilder_Tests/PangraphBuilderFromDAG_Tests/build_pangraph/"\
                "test_10_parallel_blocks_1st_and_2nd_merge_into_3rd.maf"

        expected_nodes = [
            Node(node_id=0, base=n.code('G'), aligned_to=1),
            Node(node_id=1, base=n.code('T'), aligned_to=0),
            Node(node_id=2, base=n.code('T'), aligned_to=None),
            Node(node_id=3, base=n.code('A'), aligned_to=None),
            Node(node_id=4, base=n.code('C'), aligned_to=5),
            Node(node_id=5, base=n.code('G'), aligned_to=4),
            Node(node_id=6, base=n.code('C'), aligned_to=None),

            Node(node_id=7, base=n.code('A'), aligned_to=None),
            Node(node_id=8, base=n.code('C'), aligned_to=None),
            Node(node_id=9, base=n.code('T'), aligned_to=None),
            Node(node_id=10, base=n.code('G'), aligned_to=None),
            Node(node_id=11, base=n.code('G'), aligned_to=None),

            Node(node_id=12, base=n.code('C'), aligned_to=13),
            Node(node_id=13, base=n.code('G'), aligned_to=12),
            Node(node_id=14, base=n.code('C'), aligned_to=15),
            Node(node_id=15, base=n.code('G'), aligned_to=16),
            Node(node_id=16, base=n.code('T'), aligned_to=14),
            Node(node_id=17, base=n.code('A'), aligned_to=18),
            Node(node_id=18, base=n.code('T'), aligned_to=17),
            Node(node_id=19, base=n.code('A'), aligned_to=20),
            Node(node_id=20, base=n.code('C'), aligned_to=19),
            Node(node_id=21, base=n.code('C'), aligned_to=22),
            Node(node_id=22, base=n.code('G'), aligned_to=21)
        ]

        expected_paths = {
            "seq0": [[7, 8, 9, 10, 11, 12, 15, 18, 19, 21]],
            "seq1": [[7, 8, 9, 10, 11, 12, 15, 18, 19, 21]],
            "seq2": [[0, 2, 3, 4, 6, 13, 16, 17, 20, 21]],
            "seq3": [[1, 2, 3, 5, 6, 13, 14, 17, 20, 22]]
        }

        expected_pangraph = PangraphBuilderTests.setup_pangraph(expected_nodes, expected_paths)
        actual_pangraph = PangraphBuilderTests.setup_pangraph_from_maf_firstly_converted_to_dag(maf_path,
                                                                                                self.seq_metadata,
                                                                                                self.fasta_source)
        self.compare_pangraphs(actual_pangraph=actual_pangraph, expected_pangraph=expected_pangraph)


if __name__ == '__main__':
    unittest.main()
