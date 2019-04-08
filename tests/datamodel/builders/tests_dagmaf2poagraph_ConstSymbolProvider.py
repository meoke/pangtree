import unittest
from pathlib import Path

from ...context import pNode
from ...context import pSeq
from ...context import Maf, MetadataCSV
from ...context import pPoagraph, MissingSymbol, ConstSymbolProvider
from ...context import pathtools


def nid(x): return pNode.NodeID(x)


def bid(x): return pNode.BlockID(x)


class DAGMaf2PoagraphConstSymbolProviderTests(unittest.TestCase):
    def setUp(self):
        metadata_path = Path("tests/datamodel/seq_metadata.csv")
        self.metadatacsv = MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = 'tests/datamodel/builders/maf_files_with_gaps/'
        self.missing_n = MissingSymbol()

    def test_1_missing_sequence_start(self):
        maf_path = Path(self.maf_files_dir + "test_1_missing_sequence_start.maf")
        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base(self.missing_n.value), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(1), base=pNode.Base(self.missing_n.value), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base(self.missing_n.value), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(3), base=pNode.Base('A'), aligned_to=nid(4)),
            pNode.Node(node_id=nid(4), base=pNode.Base('G'), aligned_to=nid(3)),
            pNode.Node(node_id=nid(5), base=pNode.Base('G'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(6), base=pNode.Base('G'), aligned_to=nid(7)),
            pNode.Node(node_id=nid(7), base=pNode.Base('T'), aligned_to=nid(6)),
            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(9), base=pNode.Base('A'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(10), base=pNode.Base('G'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(11), base=pNode.Base('T'), aligned_to=None, block_id=bid(0))
        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 5, 6, 11])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [4, 5, 7, 8, 9, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_2_missing_sequence_end(self):
        maf_path = Path(self.maf_files_dir + "test_2_missing_sequence_end.maf")

        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=nid(1)),
            pNode.Node(node_id=nid(1), base=pNode.Base('G'), aligned_to=nid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=nid(3)),
            pNode.Node(node_id=nid(3), base=pNode.Base('G'), aligned_to=nid(2)),
            pNode.Node(node_id=nid(4), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(5), base=pNode.Base('A'), aligned_to=nid(6)),
            pNode.Node(node_id=nid(6), base=pNode.Base('C'), aligned_to=nid(5)),

            pNode.Node(node_id=nid(7), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(8), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(10), base=pNode.Base('T'), aligned_to=None),

            pNode.Node(node_id=nid(11), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(12), base=pNode.Base(self.missing_n.value), aligned_to=None),
        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 4, 5, 8, 9, 10])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [1, 3, 4, 6, 7, 11, 12])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_3_missing_two_sequences_middle(self):
        maf_path = Path(self.maf_files_dir + "test_3_missing_two_sequences_middle.maf")

        expected_nodes = [
            # block 0
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=nid(1)),
            pNode.Node(node_id=nid(1), base=pNode.Base('G'), aligned_to=nid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=None),

            # missing seq1
            pNode.Node(node_id=nid(3), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base(self.missing_n.value), aligned_to=None),

            # missing seq2
            pNode.Node(node_id=nid(5), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(6), base=pNode.Base(self.missing_n.value), aligned_to=None),

            # block 1
            pNode.Node(node_id=nid(7), base=pNode.Base('C'), aligned_to=nid(8)),
            pNode.Node(node_id=nid(8), base=pNode.Base('G'), aligned_to=nid(7)),
            pNode.Node(node_id=nid(9), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(10), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(11), base=pNode.Base('T'), aligned_to=None)
        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 3, 4, 8, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [1, 5, 6, 7, 9, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_4_missing_one_sequence_middle(self):
        maf_path = Path(self.maf_files_dir + "test_4_missing_one_sequence_middle.maf")

        expected_nodes = [
            # block 0
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=nid(1)),
            pNode.Node(node_id=nid(1), base=pNode.Base('G'), aligned_to=nid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('A'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(5), base=pNode.Base('G'), aligned_to=nid(4)),

            # missing se2
            pNode.Node(node_id=nid(6), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(7), base=pNode.Base(self.missing_n.value), aligned_to=None),

            pNode.Node(node_id=nid(8), base=pNode.Base('A'), aligned_to=nid(9)),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=nid(8)),
            pNode.Node(node_id=nid(10), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(11), base=pNode.Base('T'), aligned_to=None),

        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 3, 4, 9, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [1, 5, 6, 7, 8, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_5_missing_one_reverted_sequence_middle_1_1(self):
        maf_path = Path(self.maf_files_dir + "test_5_missing_one_reverted_sequence_middle_1_1.maf")

        expected_nodes = [
            # block 0
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=nid(1)),
            pNode.Node(node_id=nid(1), base=pNode.Base('G'), aligned_to=nid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(4), base=pNode.Base('A'), aligned_to=nid(5)),
            pNode.Node(node_id=nid(5), base=pNode.Base('G'), aligned_to=nid(4)),

            # missing seq2, on edge (1,1)
            pNode.Node(node_id=nid(6), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(7), base=pNode.Base(self.missing_n.value), aligned_to=None),

            # block 1
            pNode.Node(node_id=nid(8), base=pNode.Base('A'), aligned_to=nid(9)),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=nid(8)),
            pNode.Node(node_id=nid(10), base=pNode.Base('G'), aligned_to=None),
            pNode.Node(node_id=nid(11), base=pNode.Base('T'), aligned_to=None)
        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 3, 4, 9, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [1, 5, 6, 7])]),
                               pSeq.SequencePath([*map(nid, [8, 10, 11])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_6_missing_one_reverted_sequence_middle_minus1_1(self):
        maf_path = Path(self.maf_files_dir + "test_6_missing_one_reverted_sequence_middle_minus1_1.maf")


        expected_nodes = [
            # block 1 because it is first in DAG and reverted
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=nid(3)),
            pNode.Node(node_id=nid(3), base=pNode.Base('T'), aligned_to=nid(2)),

            # missing seq2, on edge (-1,1)
            pNode.Node(node_id=nid(4), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(5), base=pNode.Base(self.missing_n.value), aligned_to=None),

            pNode.Node(node_id=nid(6), base=pNode.Base('A'), aligned_to=nid(7)),
            pNode.Node(node_id=nid(7), base=pNode.Base('C'), aligned_to=nid(6)),
            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(9), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(10), base=pNode.Base('A'), aligned_to=nid(11)),
            pNode.Node(node_id=nid(11), base=pNode.Base('C'), aligned_to=nid(10)),

        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2])]),
                               pSeq.SequencePath([*map(nid, [6, 8, 9, 10])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 3, 4, 5, 7, 11])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_7_missing_one_reverted_sequence_middle_minus1_minus1(self):
        maf_path = Path(self.maf_files_dir + "test_7_missing_one_reverted_sequence_middle_minus1_minus1.maf")

        expected_nodes = [
            # block 0
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None),
            pNode.Node(node_id=nid(1), base=pNode.Base('C'), aligned_to=None),
            pNode.Node(node_id=nid(2), base=pNode.Base('T'), aligned_to=None),
            pNode.Node(node_id=nid(3), base=pNode.Base('A'), aligned_to=None),

            # missing seq2
            pNode.Node(node_id=nid(4), base=pNode.Base(self.missing_n.value), aligned_to=None),
            pNode.Node(node_id=nid(5), base=pNode.Base(self.missing_n.value), aligned_to=None),

            # block 1
            pNode.Node(node_id=nid(6), base=pNode.Base('A'), aligned_to=nid(7)),
            pNode.Node(node_id=nid(7), base=pNode.Base('G'), aligned_to=nid(6)),
            pNode.Node(node_id=nid(8), base=pNode.Base('C'), aligned_to=nid(9)),
            pNode.Node(node_id=nid(9), base=pNode.Base('G'), aligned_to=nid(8)),
            pNode.Node(node_id=nid(10), base=pNode.Base('C'), aligned_to=nid(11)),
            pNode.Node(node_id=nid(11), base=pNode.Base('T'), aligned_to=nid(10)),
        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 2, 3, 7, 9, 11])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 1, 4, 5, 6, 8, 10])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph, _ = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbolProvider(self.missing_n),
            self.metadatacsv)

        self.assertEqual(expected_poagraph, actual_poagraph)

if __name__ == '__main__':
    unittest.main()
