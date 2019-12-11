import unittest
from pathlib import Path

from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.builders import po2poagraph
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.tools import pathtools


def nid(x): return graph.NodeID(x)


def bid(x): return graph.Base(x)


class Po2poagraphTests(unittest.TestCase):

    def setUp(self):
        metadata_path = Path(__file__).parent.joinpath("../seq_metadata.csv").resolve()
        self.metadatacsv = msa.MetadataCSV(pathtools.get_file_content_stringio(metadata_path), metadata_path)
        self.po_files_dir = Path(__file__).parent.joinpath("po_files").resolve()
        
    def test_1_typical_poagraph(self):
        po_path = self.po_files_dir.joinpath("test_1.po")

        expected_nodes = [graph.Node(node_id=nid(0), base=bid('A'), aligned_to=nid(1)),
                          graph.Node(node_id=nid(1), base=bid('G'), aligned_to=nid(0)),
                          graph.Node(node_id=nid(2), base=bid('C'), aligned_to=nid(3)),
                          graph.Node(node_id=nid(3), base=bid('G'), aligned_to=nid(2)),
                          graph.Node(node_id=nid(4), base=bid('A'), aligned_to=nid(5)),
                          graph.Node(node_id=nid(5), base=bid('T'), aligned_to=nid(4)),
                          graph.Node(node_id=nid(6), base=bid('G'), aligned_to=None),
                          graph.Node(node_id=nid(7), base=bid('G'), aligned_to=None),
                          graph.Node(node_id=nid(8), base=bid('A'), aligned_to=nid(9)),
                          graph.Node(node_id=nid(9), base=bid('C'), aligned_to=nid(10)),
                          graph.Node(node_id=nid(10), base=bid('G'), aligned_to=nid(11)),
                          graph.Node(node_id=nid(11), base=bid('T'), aligned_to=nid(8)),
                          graph.Node(node_id=nid(12), base=bid('A'), aligned_to=nid(13)),
                          graph.Node(node_id=nid(13), base=bid('C'), aligned_to=nid(12)),
                          graph.Node(node_id=nid(14), base=bid('T'), aligned_to=None),
                          graph.Node(node_id=nid(15), base=bid('A'), aligned_to=nid(16)),
                          graph.Node(node_id=nid(16), base=bid('C'), aligned_to=nid(17)),
                          graph.Node(node_id=nid(17), base=bid('G'), aligned_to=nid(15))
                          ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0, 2, 4, 6, 7, 8, 12, 14, 16])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [1, 2, 5, 6, 7, 9])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [3, 4, 6, 7, 10, 12, 14, 17])])],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [11, 13, 14, 15])])],
                               graph.SequenceMetadata({'group': '2'}))
        }

        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        nodes, sequences = po2poagraph.get_poagraph(msa.Po(pathtools.get_file_content_stringio(po_path), po_path),
                                                    self.metadatacsv)
        actual_poagraph = graph.Poagraph(nodes, sequences)
        self.assertEqual(expected_poagraph, actual_poagraph)

    def test_2_consensuses_and_empty_sequences(self):
        po_path = self.po_files_dir.joinpath("test_2.po")

        expected_nodes = [graph.Node(node_id=nid(0), base=bid('C'), aligned_to=nid(1)),
                          graph.Node(node_id=nid(1), base=bid('T'), aligned_to=nid(0)),
                          graph.Node(node_id=nid(2), base=bid('A'), aligned_to=nid(3)),
                          graph.Node(node_id=nid(3), base=bid('G'), aligned_to=nid(2)),
                          graph.Node(node_id=nid(4), base=bid('C'), aligned_to=None),
                          graph.Node(node_id=nid(5), base=bid('T'), aligned_to=None),
                          graph.Node(node_id=nid(6), base=bid('A'), aligned_to=nid(7)),
                          graph.Node(node_id=nid(7), base=bid('T'), aligned_to=nid(6)),
                          graph.Node(node_id=nid(8), base=bid('G'), aligned_to=None)
                          ]

        expected_sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [0, 3, 4, 5, 6, 8])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [1, 2, 4, 5, 7, 8])])],
                               graph.SequenceMetadata({'group': '1'})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [],
                               graph.SequenceMetadata({'group': '2'})),
            msa.SequenceID('CONSENS0'):
                graph.Sequence(msa.SequenceID('CONSENS0'),
                               [graph.SeqPath([*map(nid, [0, 3, 4, 5, 7, 8])])],
                               graph.SequenceMetadata({})),
            msa.SequenceID('CONSENS1'):
                graph.Sequence(msa.SequenceID('CONSENS1'),
                               [graph.SeqPath([*map(nid, [1, 2, 4, 5, 6, 8])])],
                               graph.SequenceMetadata({}))
        }

        expected_poagraph = graph.Poagraph(expected_nodes, expected_sequences)
        nodes, sequences = po2poagraph.get_poagraph(msa.Po(pathtools.get_file_content_stringio(po_path), po_path),
                                                    self.metadatacsv)
        actual_poagraph = graph.Poagraph(nodes, sequences)
        self.assertEqual(expected_poagraph, actual_poagraph)


if __name__ == '__main__':
    unittest.main()
