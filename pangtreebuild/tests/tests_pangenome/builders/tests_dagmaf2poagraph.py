import unittest
from pathlib import Path

from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.pangenome import builder
from pangtreebuild.pangenome.parameters import missings
from pangtreebuild.tools import pathtools


def nid(x): return graph.NodeID(x)


def bid(x): return graph.BlockID(x)


class DAGMaf2PoagraphFakeFastaProviderTests(unittest.TestCase):

    def setUp(self):
        metadata_path = Path(__file__).parent.joinpath("../seq_metadata.csv").resolve()
        self.metadatacsv = msa.MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = Path(__file__).parent.joinpath("maf_files_with_cycles_or_reversion").resolve()
        self.fasta_provider = missings.ConstBaseProvider(missings.MissingBase())

    def test_00_simple(self):
        maf_path = self.maf_files_dir.joinpath("test_0_simple.maf")
        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(3), base=graph.Base('G'), aligned_to=None, block_id=bid(0)),

            graph.Node(node_id=nid(4), base=graph.Base('A'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(5), base=graph.Base('C'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(6), base=graph.Base('T'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(7), base=graph.Base('G'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(8), base=graph.Base('A'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(9), base=graph.Base('A'), aligned_to=None, block_id=bid(1))
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 4, 5, 6, 7, 8])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_01_reversed_seq_in_one_block(self):
        maf_path = self.maf_files_dir.joinpath("test_1_reversed_seq_in_one_block.maf")
        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=nid(1), block_id=bid(0)),
            graph.Node(node_id=nid(1), base=graph.Base('G'), aligned_to=nid(0), block_id=bid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=nid(3), block_id=bid(0)),
            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=nid(2), block_id=bid(0)),

            graph.Node(node_id=nid(4), base=graph.Base('G'), aligned_to=nid(5), block_id=bid(1)),
            graph.Node(node_id=nid(5), base=graph.Base('T'), aligned_to=nid(4), block_id=bid(1)),
            graph.Node(node_id=nid(6), base=graph.Base('A'), aligned_to=nid(7), block_id=bid(1)),
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=nid(6), block_id=bid(1)),
            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=None, block_id=bid(1)),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0, 2, 4, 7, 8])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [1, 3])]),
                                graph.SeqPath([*map(nid, [5, 6])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_02_seq_starts_in_second_block(self):
        maf_path = self.maf_files_dir.joinpath(
                        "test_2_seq_starts_in_second_block.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('C'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(1), base=graph.Base('T'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('G'), aligned_to=None, block_id=bid(0)),

            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=None, block_id=bid(1)),

            graph.Node(node_id=nid(4), base=graph.Base('G'), aligned_to=nid(5), block_id=bid(2)),
            graph.Node(node_id=nid(5), base=graph.Base('T'), aligned_to=nid(4), block_id=bid(2)),
            graph.Node(node_id=nid(6), base=graph.Base('A'), aligned_to=None, block_id=bid(2)),
            graph.Node(node_id=nid(7), base=graph.Base('A'), aligned_to=None, block_id=bid(2)),
            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=None, block_id=bid(2)),

        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [1, 2, 3])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 5, 7])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [3, 4, 6, 8])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_03_edge_not_from_last_node_in_block(self):
        maf_path = self.maf_files_dir.joinpath("test_3_edge_not_from_last_node_in_block.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(5), base=graph.Base('A'), aligned_to=nid(6)),
            graph.Node(node_id=nid(6), base=graph.Base('T'), aligned_to=nid(5)),
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=None),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 3, 4, 6, 7])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 5])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_04_single_block_no_nucleotides(self):
        maf_path = self.maf_files_dir.joinpath(
                        "test_4_single_block_no_nucleotides.maf")

        expected_nodes = []

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_05_single_block_single_nucletodide(self):
        maf_path = self.maf_files_dir.joinpath(
                        "test_5_single_block_single_nucletodide.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None, block_id=bid(0))
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [0])])],
                               graph.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_06_1st_block_separates_into_2_branches_which_connect_in_3rd_block(self):
        maf_path = self.maf_files_dir.joinpath(
                        "test_6_1st_block_separates_into_2_branches_which_connect_in_3rd_block.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=nid(1)),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=nid(2)),
            graph.Node(node_id=nid(2), base=graph.Base('G'), aligned_to=nid(0)),
            graph.Node(node_id=nid(3), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('A'), aligned_to=nid(5)),
            graph.Node(node_id=nid(5), base=graph.Base('T'), aligned_to=nid(4)),

            graph.Node(node_id=nid(6), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(7), base=graph.Base('G'), aligned_to=None),

            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=nid(9)),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=nid(10)),
            graph.Node(node_id=nid(10), base=graph.Base('T'), aligned_to=nid(8)),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0, 3, 4, 8])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [1, 3, 5, 6, 7, 9])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [2, 3, 5, 10])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_07_inactive_edges_due_to_reversed_seqs(self):
        maf_path = self.maf_files_dir.joinpath("test_7_inactive_edges_due_to_reversed_seqs.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('A'), aligned_to=nid(5)),
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=nid(4)),
            graph.Node(node_id=nid(6), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(7), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(8), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(10), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(11), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(12), base=graph.Base('A'), aligned_to=nid(13)),
            graph.Node(node_id=nid(13), base=graph.Base('C'), aligned_to=nid(12)),
            graph.Node(node_id=nid(14), base=graph.Base('A'), aligned_to=nid(15)),
            graph.Node(node_id=nid(15), base=graph.Base('T'), aligned_to=nid(14)),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 4, 6, 8, 9, 10, 11, 12, 14])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 1, 2])]),
                                graph.SeqPath([*map(nid, [3, 5])]),
                                graph.SeqPath([*map(nid, [6, 7, 8, 9, 10, 11, 13, 15])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [0, 2])]),
                                graph.SeqPath([*map(nid, [3, 5])]),
                                graph.SeqPath([*map(nid, [6, 7, 8, 9, 11, 12, 15])])],
                               graph.SequenceMetadata({'group': '2'})),
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_08_reversed_block(self):
        maf_path = self.maf_files_dir.joinpath("test_8_reversed_block.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(1), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None),
            # next block is reversed because it was converted to dag
            graph.Node(node_id=nid(3), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(5), base=graph.Base('A'), aligned_to=nid(6)),
            graph.Node(node_id=nid(6), base=graph.Base('G'), aligned_to=nid(5)),
            graph.Node(node_id=nid(7), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(8), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(9), base=graph.Base('T'), aligned_to=None),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 3, 4, 5, 7, 8, 9])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 4, 6, 7, 8, 9])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 4, 6, 7, 9])])],
                               graph.SequenceMetadata({'group': '2'})),
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_09_inactive_edges_but_all_strands_plus(self):
        maf_path = self.maf_files_dir.joinpath("test_9_inactive_edges_but_all_strands_plus.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('G'), aligned_to=None),

            graph.Node(node_id=nid(5), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(6), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(7), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(8), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=None),

            graph.Node(node_id=nid(10), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(11), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(12), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(13), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(14), base=graph.Base('G'), aligned_to=None),

            graph.Node(node_id=nid(15), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(16), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(17), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(18), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(19), base=graph.Base('G'), aligned_to=None),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 4, 10, 11, 12, 13, 14])]),
                                graph.SeqPath([*map(nid, [5, 6, 7, 8, 9, 15, 16, 17, 18, 19])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                          13, 14, 15, 16, 17, 18, 19])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_10_parallel_blocks_1st_and_2nd_merge_into_3rd(self):
        maf_path = self.maf_files_dir.joinpath("test_10_parallel_blocks_1st_and_2nd_merge_into_3rd.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('G'), aligned_to=nid(1)),
            graph.Node(node_id=nid(1), base=graph.Base('T'), aligned_to=nid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('C'), aligned_to=nid(5)),
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=nid(4)),
            graph.Node(node_id=nid(6), base=graph.Base('C'), aligned_to=None),

            graph.Node(node_id=nid(7), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(9), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(10), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(11), base=graph.Base('G'), aligned_to=None),

            graph.Node(node_id=nid(12), base=graph.Base('C'), aligned_to=nid(13)),
            graph.Node(node_id=nid(13), base=graph.Base('G'), aligned_to=nid(12)),
            graph.Node(node_id=nid(14), base=graph.Base('C'), aligned_to=nid(15)),
            graph.Node(node_id=nid(15), base=graph.Base('G'), aligned_to=nid(16)),
            graph.Node(node_id=nid(16), base=graph.Base('T'), aligned_to=nid(14)),
            graph.Node(node_id=nid(17), base=graph.Base('A'), aligned_to=nid(18)),
            graph.Node(node_id=nid(18), base=graph.Base('T'), aligned_to=nid(17)),
            graph.Node(node_id=nid(19), base=graph.Base('A'), aligned_to=nid(20)),
            graph.Node(node_id=nid(20), base=graph.Base('C'), aligned_to=nid(19)),
            graph.Node(node_id=nid(21), base=graph.Base('C'), aligned_to=nid(22)),
            graph.Node(node_id=nid(22), base=graph.Base('G'), aligned_to=nid(21)),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [7, 8, 9, 10, 11, 12, 15, 18, 19, 21])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [7, 8, 9, 10, 11, 12, 15, 18, 19, 21])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 2, 3, 4, 6, 13, 16, 17, 20, 21])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [1, 2, 3, 5, 6, 13, 14, 17, 20, 22])])],
                               graph.SequenceMetadata({'group': '2'})),
        }
        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)


if __name__ == '__main__':
    unittest.main()
