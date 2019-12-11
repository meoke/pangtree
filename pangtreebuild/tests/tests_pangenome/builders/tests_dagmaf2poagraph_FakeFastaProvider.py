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

    class FakeFastaProvider(missings.FastaProvider):
        def __init__(self):
            self.sources = {
                msa.SequenceID("seq0"): "",
                msa.SequenceID("seq1"): "ACTAGGT",
                msa.SequenceID("seq2"): "GGTCAGT",
                msa.SequenceID("seq3"): "",
                msa.SequenceID("seq4"): ""}

        def get_base(self, sequence_id: msa.SequenceID, i: int) -> graph.Base:
            try:
                return graph.Base(self.sources[sequence_id][i])
            except KeyError:
                raise Exception("No record found with given id_!")

    def setUp(self):
        metadata_path = Path(__file__).parent.joinpath("../seq_metadata.csv").resolve()
        self.metadatacsv = msa.MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = Path(__file__).parent.joinpath("maf_files_with_gaps").resolve()
        self.fasta_provider = DAGMaf2PoagraphFakeFastaProviderTests.FakeFastaProvider()

    def test_1_missing_sequence_start(self):
        maf_path = self.maf_files_dir.joinpath("test_1_missing_sequence_start.maf")
        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(3), base=graph.Base('A'), aligned_to=nid(4)),
            graph.Node(node_id=nid(4), base=graph.Base('G'), aligned_to=nid(3)),
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(6), base=graph.Base('G'), aligned_to=nid(7)),
            graph.Node(node_id=nid(7), base=graph.Base('T'), aligned_to=nid(6)),
            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(9), base=graph.Base('A'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(10), base=graph.Base('G'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(11), base=graph.Base('T'), aligned_to=None, block_id=bid(0))
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 5, 6, 11])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [4, 5, 7, 8, 9, 10, 11])])],
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

    def test_2_missing_sequence_end(self):
        maf_path = self.maf_files_dir.joinpath("test_2_missing_sequence_end.maf")

        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=nid(1)),
            graph.Node(node_id=nid(1), base=graph.Base('G'), aligned_to=nid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=nid(3)),
            graph.Node(node_id=nid(3), base=graph.Base('G'), aligned_to=nid(2)),
            graph.Node(node_id=nid(4), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(5), base=graph.Base('A'), aligned_to=nid(6)),
            graph.Node(node_id=nid(6), base=graph.Base('C'), aligned_to=nid(5)),

            graph.Node(node_id=nid(7), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(8), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(10), base=graph.Base('T'), aligned_to=None),

            graph.Node(node_id=nid(11), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(12), base=graph.Base('T'), aligned_to=None),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 2, 4, 5, 8, 9, 10])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [1, 3, 4, 6, 7, 11, 12])])],
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

    def test_3_missing_two_sequences_middle(self):
        maf_path = self.maf_files_dir.joinpath("test_3_missing_two_sequences_middle.maf")

        expected_nodes = [
            # block 0
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=nid(1)),
            graph.Node(node_id=nid(1), base=graph.Base('G'), aligned_to=nid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=None),

            # missing seq1
            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('A'), aligned_to=None),

            # missing seq2
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(6), base=graph.Base('T'), aligned_to=None),

            # block 1
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=nid(8)),
            graph.Node(node_id=nid(8), base=graph.Base('G'), aligned_to=nid(7)),
            graph.Node(node_id=nid(9), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(10), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(11), base=graph.Base('T'), aligned_to=None)
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 2, 3, 4, 8, 10, 11])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [1, 5, 6, 7, 9, 10, 11])])],
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

    def test_4_missing_one_sequence_middle(self):
        maf_path = self.maf_files_dir.joinpath("test_4_missing_one_sequence_middle.maf")

        expected_nodes = [
            # block 0
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=nid(1)),
            graph.Node(node_id=nid(1), base=graph.Base('G'), aligned_to=nid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('A'), aligned_to=nid(5)),
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=nid(4)),

            # missing se2
            graph.Node(node_id=nid(6), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=None),

            graph.Node(node_id=nid(8), base=graph.Base('A'), aligned_to=nid(9)),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=nid(8)),
            graph.Node(node_id=nid(10), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(11), base=graph.Base('T'), aligned_to=None),

        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 2, 3, 4, 9, 10, 11])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [1, 5, 6, 7, 8, 10, 11])])],
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

    def test_5_missing_one_reverted_sequence_middle_1_1(self):
        maf_path = self.maf_files_dir.joinpath("test_5_missing_one_reverted_sequence_middle_1_1.maf")

        expected_nodes = [
            # block 0
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=nid(1)),
            graph.Node(node_id=nid(1), base=graph.Base('G'), aligned_to=nid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(4), base=graph.Base('A'), aligned_to=nid(5)),
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=nid(4)),

            # missing seq2, on edge (1,1)
            graph.Node(node_id=nid(6), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=None),

            # block 1
            graph.Node(node_id=nid(8), base=graph.Base('A'), aligned_to=nid(9)),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=nid(8)),
            graph.Node(node_id=nid(10), base=graph.Base('G'), aligned_to=None),
            graph.Node(node_id=nid(11), base=graph.Base('T'), aligned_to=None)
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 2, 3, 4, 9, 10, 11])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [1, 5, 6, 7])]),
                                graph.SeqPath([*map(nid, [8, 10, 11])])],
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

    def test_6_missing_one_reverted_sequence_middle_minus1_1(self):
        maf_path = self.maf_files_dir.joinpath("test_6_missing_one_reverted_sequence_middle_minus1_1.maf")

        expected_nodes = [
            # block 1 because it is first in DAG and reverted
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=nid(3)),
            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=nid(2)),

            # missing seq2, on edge (-1,1)
            graph.Node(node_id=nid(4), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(5), base=graph.Base('C'), aligned_to=None),

            graph.Node(node_id=nid(6), base=graph.Base('A'), aligned_to=nid(7)),
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=nid(6)),
            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(9), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(10), base=graph.Base('A'), aligned_to=nid(11)),
            graph.Node(node_id=nid(11), base=graph.Base('C'), aligned_to=nid(10)),

        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 2])]),
                                graph.SeqPath([*map(nid, [6, 8, 9, 10])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 1, 3, 4, 5, 7, 11])])],
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

    def test_7_missing_one_reverted_sequence_middle_minus1_minus1(self):
        maf_path = self.maf_files_dir.joinpath("test_7_missing_one_reverted_sequence_middle_minus1_minus1.maf")

        expected_nodes = [
            # block 0
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None),
            graph.Node(node_id=nid(1), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(2), base=graph.Base('T'), aligned_to=None),
            graph.Node(node_id=nid(3), base=graph.Base('A'), aligned_to=None),

            # missing seq2
            graph.Node(node_id=nid(4), base=graph.Base('C'), aligned_to=None),
            graph.Node(node_id=nid(5), base=graph.Base('A'), aligned_to=None),

            # block 1
            graph.Node(node_id=nid(6), base=graph.Base('A'), aligned_to=nid(7)),
            graph.Node(node_id=nid(7), base=graph.Base('G'), aligned_to=nid(6)),
            graph.Node(node_id=nid(8), base=graph.Base('C'), aligned_to=nid(9)),
            graph.Node(node_id=nid(9), base=graph.Base('G'), aligned_to=nid(8)),
            graph.Node(node_id=nid(10), base=graph.Base('C'), aligned_to=nid(11)),
            graph.Node(node_id=nid(11), base=graph.Base('T'), aligned_to=nid(10)),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [0, 1, 2, 3, 7, 9, 11])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 1, 4, 5, 6, 8, 10])])],
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


if __name__ == '__main__':
    unittest.main()
