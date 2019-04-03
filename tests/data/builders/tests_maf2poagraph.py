import unittest
from pathlib import Path

from ...context import pNode
from ...context import pSeq
from ...context import Maf, MetadataCSV
from ...context import maf2poagraph, pathtools


def nid(x): return pNode.NodeID(x)


def bid(x): return pNode.BlockID(x)


class Maf2poagraphTests(unittest.TestCase):

    def setUp(self):
        metadata_path = Path("tests/data/seq_metadata.csv")
        self.metadatacsv = MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = 'tests/data/builders/maf_files/'

    def test_1_messy_sequences(self):
        maf_path = Path(self.maf_files_dir + "test_1_messy_sequences.maf")
        expected_nodes = [
            pNode.Node(node_id=nid(0), base=pNode.Base('A'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(1), base=pNode.Base('A'), aligned_to=nid(2), block_id=bid(0)),
            pNode.Node(node_id=nid(2), base=pNode.Base('C'), aligned_to=nid(1), block_id=bid(0)),
            pNode.Node(node_id=nid(3), base=pNode.Base('T'), aligned_to=None, block_id=bid(0)),
            pNode.Node(node_id=nid(4), base=pNode.Base('C'), aligned_to=nid(5), block_id=bid(0)),
            pNode.Node(node_id=nid(5), base=pNode.Base('G'), aligned_to=nid(4), block_id=bid(0)),
            pNode.Node(node_id=nid(6), base=pNode.Base('A'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(7), base=pNode.Base('C'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(8), base=pNode.Base('G'), aligned_to=None, block_id=bid(1)),
            pNode.Node(node_id=nid(9), base=pNode.Base('C'), aligned_to=nid(10), block_id=bid(2)),
            pNode.Node(node_id=nid(10), base=pNode.Base('G'), aligned_to=nid(9), block_id=bid(2)),
            pNode.Node(node_id=nid(11), base=pNode.Base('T'), aligned_to=None, block_id=bid(2)),
            pNode.Node(node_id=nid(12), base=pNode.Base('C'), aligned_to=None, block_id=bid(2)),
            pNode.Node(node_id=nid(13), base=pNode.Base('C'), aligned_to=None, block_id=bid(2)),
            pNode.Node(node_id=nid(14), base=pNode.Base('A'), aligned_to=None, block_id=bid(2)),
        ]

        expected_sequences = {
            pSeq.SequenceID('seq0'):
                pSeq.Sequence(pSeq.SequenceID('seq0'),
                              [pSeq.SequencePath([*map(nid, [1, 3, 4, 6, 8, 9, 11, 12])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq1'):
                pSeq.Sequence(pSeq.SequenceID('seq1'),
                              [pSeq.SequencePath([*map(nid, [2, 3, 4, 10, 11, 12, 13, 14])])],
                              pSeq.SequenceMetadata({'group': '1'})),
            pSeq.SequenceID('seq2'):
                pSeq.Sequence(pSeq.SequenceID('seq2'),
                              [pSeq.SequencePath([*map(nid, [0, 2, 5, 6, 7, 10, 11, 12, 14])])],
                              pSeq.SequenceMetadata({'group': '2'})),
            pSeq.SequenceID('seq3'):
                pSeq.Sequence(pSeq.SequenceID('seq3'),
                              [],
                              pSeq.SequenceMetadata({'group': '2'}))
        }
        actual_nodes, actual_sequences = maf2poagraph.get_poagraph(Maf(pathtools.get_file_content_stringio(maf_path),
                                                                       maf_path),
                                                                   self.metadatacsv)

        self.assertEqual(expected_nodes, actual_nodes)
        self.assertEqual(expected_sequences, actual_sequences)


if __name__ == '__main__':
    unittest.main()
