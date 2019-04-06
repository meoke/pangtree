import unittest
from pathlib import Path

from ...context import pNode
from ...context import pSeq
from ...context import Maf, MetadataCSV
from ...context import pPoagraph, ConstSymbolProvider, MissingSymbol
from ...context import pathtools


def nid(x): return pNode.NodeID(x)


def bid(x): return pNode.BlockID(x)


class DAGMaf2PoagraphFakeFastaProviderTests(unittest.TestCase):

    def setUp(self):
        metadata_path = Path("tests/datamodel/seq_metadata.csv")
        self.metadatacsv = MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = 'tests/datamodel/builders/maf_files_with_cycles_or_reversion/'
        self.fasta_provider = ConstSymbolProvider(MissingSymbol())

    def test_00_simple(self):
        maf_path = Path(self.maf_files_dir + "test_0_simple.maf")
        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(3), base=pNode.Base('G'), aligned_to=None, block_id=bid(0)),

            pNode.Node(node_id=nid(4), base=pNode.Base('A'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(5), base=pNode.Base('C'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(6), base=pNode.Base('T'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(7), base=pNode.Base('G'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(8), base=pNode.Base('A'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(9), base=pNode.Base('A'), aligned_to=None, block_id=bid(1))
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4, 5, 6, 7, 8])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'})),
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_01_reversed_seq_in_one_block(self):
        maf_path = Path(self.maf_files_dir + "test_1_reversed_seq_in_one_block.maf")
        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=nid(1), block_id=bid(0)),
            pNode.Node(node_id=nid(1), base=pNode.Base('G'), aligned_to=nid(0), block_id=bid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=nid(3), block_id=bid(0)),
            pNode.Node(node_id=nid(3), base=pNode.Base('T'), aligned_to=nid(2), block_id=bid(0)),

            pNode.Node(node_id=nid(4), base=pNode.Base('G'), aligned_to=nid(5), block_id=bid(1)),
            pNode.Node(node_id=nid(5), base=pNode.Base('T'), aligned_to=nid(4), block_id=bid(1)),
            pNode.Node(node_id=nid(6), base=pNode.Base('A'), aligned_to=nid(7), block_id=bid(1)),
            pNode.Node(node_id=nid(7), base=pNode.Base('C'), aligned_to=nid(6), block_id=bid(1)),
            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=None, block_id=bid(1)),
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 4, 7, 8])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [1, 3])]),
                               pSeq.SequencePath([*map(nid, [5, 6])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'})),
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_02_seq_starts_in_second_block(self):
        maf_path = Path(self.maf_files_dir + "test_2_seq_starts_in_second_block.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('C'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(1), base=pNode.Base('T'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('G'), aligned_to=None, block_id=bid(0)),

            pNode.Node(node_id=nid(3), base=pNode.Base('T'), aligned_to=None, block_id=bid(1)),

            pNode.Node(node_id=nid(4), base=pNode.Base('G'), aligned_to=nid(5), block_id=bid(2)),
            pNode.Node(node_id=nid(5), base=pNode.Base('T'), aligned_to=nid(4), block_id=bid(2)),
            pNode.Node(node_id=nid(6), base=pNode.Base('A'), aligned_to=None, block_id=bid(2)),
            pNode.Node(node_id=nid(7), base=pNode.Base('A'), aligned_to=None, block_id=bid(2)),
            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=None, block_id=bid(2)),

        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [1, 2, 3])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 5, 7])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [3, 4, 6, 8])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_03_edge_not_from_last_node_in_block(self):
        maf_path = Path(self.maf_files_dir + "test_3_edge_not_from_last_node_in_block.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(5), base=pNode.Base('A'), aligned_to=nid(6)),
            pNode.Node(node_id=nid(6), base=pNode.Base('T'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(7), base=pNode.Base('C'), aligned_to=None),
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 3, 4, 6, 7])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 5])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_04_single_block_no_nucleotides(self):
        maf_path = Path(self.maf_files_dir + "test_4_single_block_no_nucleotides.maf")

        expected_nodes = []

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_05_single_block_single_nucletodide(self):
        maf_path = Path(self.maf_files_dir + "test_5_single_block_single_nucletodide.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None, block_id=bid(0))
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [0])])],
                              pSeq.SequenceMetadata({'group': '2'}))
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_06_1st_block_separates_into_2_branches_which_connect_in_3rd_block(self):
        maf_path = Path(self.maf_files_dir + "test_6_1st_block_separates_into_2_branches_which_connect_in_3rd_block.maf")


        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=nid(1)),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=nid(2)),
            pNode.Node(node_id=nid(2), base=pNode.Base('G'), aligned_to=nid(0)),
            pNode.Node(node_id=nid(3), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('A'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(5), base=pNode.Base('T'), aligned_to=nid(4)),

            pNode.Node(node_id=nid(6), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(7), base=pNode.Base('G'), aligned_to=None),

            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=nid(9)),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=nid(10)),
            pNode.Node(node_id=nid(10), base=pNode.Base('T'), aligned_to=nid(8)),

        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [0, 3, 4, 8])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [1, 3, 5, 6, 7, 9])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [2, 3, 5, 10])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_07_inactive_edges_due_to_reversed_seqs(self):
        maf_path = Path(self.maf_files_dir + "test_7_inactive_edges_due_to_reversed_seqs.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('A'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(5), base=pNode.Base('G'), aligned_to=nid(4)),
            pNode.Node(node_id=nid(6), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(7), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(8), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(10), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(11), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(12), base=pNode.Base('A'), aligned_to=nid(13)),
            pNode.Node(node_id=nid(13), base=pNode.Base('C'), aligned_to=nid(12)),
            pNode.Node(node_id=nid(14), base=pNode.Base('A'), aligned_to=nid(15)),
            pNode.Node(node_id=nid(15), base=pNode.Base('T'), aligned_to=nid(14)),
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 4, 6, 8, 9, 10, 11, 12, 14])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2])]),
                               pSeq.SequencePath([*map(nid, [3, 5])]),
                               pSeq.SequencePath([*map(nid, [6, 7, 8, 9, 10, 11, 13, 15])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [0, 2])]),
                               pSeq.SequencePath([*map(nid, [3, 5])]),
                               pSeq.SequencePath([*map(nid, [6, 7, 8, 9, 11, 12, 15])])],
                              pSeq.SequenceMetadata({'group': '2'})),
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_08_reversed_block(self):
        maf_path = Path(self.maf_files_dir + "test_8_reversed_block.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(1), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None),
            # next block is reversed because it was converted to dag
            pNode.Node(node_id=nid(3), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(5), base=pNode.Base('A'), aligned_to=nid(6)),
            pNode.Node(node_id=nid(6), base=pNode.Base('G'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(7), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(8), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(9), base=pNode.Base('T'), aligned_to=None),
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 3, 4, 5, 7, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4, 6, 7, 8, 9])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4, 6, 7, 9])])],
                              pSeq.SequenceMetadata({'group': '2'})),
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_09_inactive_edges_but_all_strands_plus(self):
        maf_path = Path(self.maf_files_dir + "test_9_inactive_edges_but_all_strands_plus.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('G'), aligned_to=None),

            pNode.Node(node_id=nid(5), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(6), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(7), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(8), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=None),

            pNode.Node(node_id=nid(10), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(11), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(12), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(13), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(14), base=pNode.Base('G'), aligned_to=None),

            pNode.Node(node_id=nid(15), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(16), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(17), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(18), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(19), base=pNode.Base('G'), aligned_to=None),
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4, 10, 11, 12, 13, 14])]),
                               pSeq.SequencePath([*map(nid, [5, 6, 7, 8, 9, 15, 16, 17, 18, 19])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'})),
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_10_parallel_blocks_1st_and_2nd_merge_into_3rd(self):
        maf_path = Path(self.maf_files_dir + "test_10_parallel_blocks_1st_and_2nd_merge_into_3rd.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('G'), aligned_to=nid(1)),
            pNode.Node(node_id=nid(1), base=pNode.Base('T'), aligned_to=nid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('C'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(5), base=pNode.Base('G'), aligned_to=nid(4)),
            pNode.Node(node_id=nid(6), base=pNode.Base('C'), aligned_to=None),

            pNode.Node(node_id=nid(7), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(9), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(10), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(11), base=pNode.Base('G'), aligned_to=None),

            pNode.Node(node_id=nid(12), base=pNode.Base('C'), aligned_to=nid(13)),
            pNode.Node(node_id=nid(13), base=pNode.Base('G'), aligned_to=nid(12)),
            pNode.Node(node_id=nid(14), base=pNode.Base('C'), aligned_to=nid(15)),
            pNode.Node(node_id=nid(15), base=pNode.Base('G'), aligned_to=nid(16)),
            pNode.Node(node_id=nid(16), base=pNode.Base('T'), aligned_to=nid(14)),
            pNode.Node(node_id=nid(17), base=pNode.Base('A'), aligned_to=nid(18)),
            pNode.Node(node_id=nid(18), base=pNode.Base('T'), aligned_to=nid(17)),
            pNode.Node(node_id=nid(19), base=pNode.Base('A'), aligned_to=nid(20)),
            pNode.Node(node_id=nid(20), base=pNode.Base('C'), aligned_to=nid(19)),
            pNode.Node(node_id=nid(21), base=pNode.Base('C'), aligned_to=nid(22)),
            pNode.Node(node_id=nid(22), base=pNode.Base('G'), aligned_to=nid(21)),
        ]

        expected_sequences = pSeq.Sequences({
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [7, 8, 9, 10, 11, 12, 15, 18, 19, 21])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [7, 8, 9, 10, 11, 12, 15, 18, 19, 21])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 3, 4, 6, 13, 16, 17, 20, 21])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [pSeq.SequencePath([*map(nid, [1, 2, 3, 5, 6, 13, 14, 17, 20, 22])])],
                              pSeq.SequenceMetadata({'group': '2'})),
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

if __name__ == '__main__':
    unittest.main()
