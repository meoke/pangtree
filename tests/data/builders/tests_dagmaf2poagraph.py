import unittest
from pathlib import Path

from ...context import pNode
from ...context import pSeq
from ...context import Maf, MetadataCSV
from ...context import pPoagraph, MissingSymbol, ConstSymbol
from ...context import pathtools


def nid(x): return pNode.NodeID(x)


def bid(x): return pNode.BlockID(x)


class DAGMaf2PoagraphTests(unittest.TestCase):
    def setUp(self):
        metadata_path = Path("tests/data/seq_metadata.csv")
        self.metadatacsv = MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.maf_files_dir = 'tests/data/fasta_providers/maf_files_with_gaps'
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

        expected_sequences = pSeq.Sequences({
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
        })
        expected_poagraph = pPoagraph.Poagraph(expected_nodes, expected_sequences)
        actual_poagraph = pPoagraph.Poagraph.build_from_dagmaf(
            Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            ConstSymbol(self.missing_n),
            self.metadatacsv)
        self.assertEqual(expected_poagraph, actual_poagraph)


if __name__ == '__main__':
    unittest.main()
