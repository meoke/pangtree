import unittest
from pathlib import Path

from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.builders import maf2poagraph
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.tools import pathtools


def nid(x): return graph.NodeID(x)


def bid(x): return graph.BlockID(x)


class Maf2poagraphTests(unittest.TestCase):

    def setUp(self):
        metadata_path = Path(__file__).parent.joinpath("../seq_metadata.csv").resolve()
        self.metadatacsv = msa.MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = Path(__file__).parent.joinpath("maf_files").resolve()

    def test_1_messy_sequences(self):
        maf_path = self.maf_files_dir.joinpath("test_1_messy_sequences.maf")
        expected_nodes = [
            graph.Node(node_id=nid(0), base=graph.Base('A'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(1), base=graph.Base('A'), aligned_to=nid(2), block_id=bid(0)),
            graph.Node(node_id=nid(2), base=graph.Base('C'), aligned_to=nid(1), block_id=bid(0)),
            graph.Node(node_id=nid(3), base=graph.Base('T'), aligned_to=None, block_id=bid(0)),
            graph.Node(node_id=nid(4), base=graph.Base('C'), aligned_to=nid(5), block_id=bid(0)),
            graph.Node(node_id=nid(5), base=graph.Base('G'), aligned_to=nid(4), block_id=bid(0)),
            graph.Node(node_id=nid(6), base=graph.Base('A'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(7), base=graph.Base('C'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(8), base=graph.Base('G'), aligned_to=None, block_id=bid(1)),
            graph.Node(node_id=nid(9), base=graph.Base('C'), aligned_to=nid(10), block_id=bid(2)),
            graph.Node(node_id=nid(10), base=graph.Base('G'), aligned_to=nid(9), block_id=bid(2)),
            graph.Node(node_id=nid(11), base=graph.Base('T'), aligned_to=None, block_id=bid(2)),
            graph.Node(node_id=nid(12), base=graph.Base('C'), aligned_to=None, block_id=bid(2)),
            graph.Node(node_id=nid(13), base=graph.Base('C'), aligned_to=None, block_id=bid(2)),
            graph.Node(node_id=nid(14), base=graph.Base('A'), aligned_to=None, block_id=bid(2)),
        ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [1, 3, 4, 6, 8, 9, 11, 12])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [2, 3, 4, 10, 11, 12, 13, 14])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [0, 2, 5, 6, 7, 10, 11, 12, 14])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'}))
        }
        actual_nodes, actual_sequences = maf2poagraph.get_poagraph(msa.Maf(pathtools.get_file_content_stringio(maf_path),
                                                                           maf_path),
                                                                   self.metadatacsv)

        self.assertEqual(expected_nodes, actual_nodes)
        self.assertEqual(expected_sequences, actual_sequences)


if __name__ == '__main__':
    unittest.main()
