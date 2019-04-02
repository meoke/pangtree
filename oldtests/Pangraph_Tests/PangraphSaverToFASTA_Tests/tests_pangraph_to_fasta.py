import unittest
from ddt import ddt

from tests.Pangraph_Tests.Pangraph_Tests import PangraphTests

from tests.context import make_base as n
from tests.context import Node, NodeID
from tests.context import MultialignmentMetadata
from tests.context import fasta
from tests.context import ConsensusesTree, ConsensusNode, ConsensusNodeID, SequenceID
from tests.context import CompatibilityToPath
from tools import pathtools


@ddt
class PangraphSaverToFASTATest_BuildPangraph(PangraphTests):

    def setUp(self):
        self.seq_metadata = MultialignmentMetadata(
            pathtools.get_file_content("Pangraph_Tests/seq_metadata.csv"))

        pangraph_nodes = [Node(node_id=0, base=n('A'), aligned_to=1),
                          Node(node_id=1, base=n('G'), aligned_to=0),
                          Node(node_id=2, base=n('C'), aligned_to=3),
                          Node(node_id=3, base=n('G'), aligned_to=2),
                          Node(node_id=4, base=n('A'), aligned_to=5),
                          Node(node_id=5, base=n('T'), aligned_to=4),
                          Node(node_id=6, base=n('G'), aligned_to=None),
                          Node(node_id=7, base=n('G'), aligned_to=None),
                          Node(node_id=8, base=n('A'), aligned_to=9),
                          Node(node_id=9, base=n('C'), aligned_to=10),
                          Node(node_id=10, base=n('G'), aligned_to=11),
                          Node(node_id=11, base=n('T'), aligned_to=8),
                          Node(node_id=12, base=n('A'), aligned_to=13),
                          Node(node_id=13, base=n('C'), aligned_to=12),
                          Node(node_id=14, base=n('T'), aligned_to=None),
                          Node(node_id=15, base=n('A'), aligned_to=16),
                          Node(node_id=16, base=n('C'), aligned_to=17),
                          Node(node_id=17, base=n('G'), aligned_to=15)]

        pangraph_paths = {
            'seq0': [[0, 2, 4, 6, 7, 8, 12, 14, 16]],
            'seq1': [],
            'seq2': [[3, 4, 6, 7, 10, 12], [14, 17]],
            'seq3': [[11], [13, 14, 15]]
        }

        self.pangraph = PangraphTests.setup_pangraph(pangraph_nodes, pangraph_paths)

    def test_1_sequences_fasta(self):
        expected_sequences_fasta_path = "Pangraph_Tests/" \
                                                "PangraphSaverToFASTA_Tests/" \
                                                "files/" \
                                                "sequences.fasta"

        actual_sequences_fasta_content = fasta.pangraph_to_fasta(self.pangraph)
        expected_sequences_fasta_content = PangraphTests.get_file_content_as_str(expected_sequences_fasta_path)
        self.assertEqual(expected_sequences_fasta_content, actual_sequences_fasta_content)


    def test_2_consensuses_tree_fasta(self):
        expected_consensuses_fasta_path = "Pangraph_Tests/" \
                                                  "PangraphSaverToFASTA_Tests/" \
                                                  "files/" \
                                                  "consensuses.fasta"
        consensuses_tree = ConsensusesTree()
        consensuses_tree.nodes = [
            # all members set
            ConsensusNode(parent_node_id=ConsensusNodeID(-1),
                          children_nodes_ids=[ConsensusNodeID(1), ConsensusNodeID(2)],
                          sequences_ids=[SequenceID('seq0'), SequenceID('seq1'), SequenceID('seq2'),
                                         SequenceID('seq3')],
                          consensus_id=ConsensusNodeID(0),
                          mincomp=CompatibilityToPath(0.5, 1),
                          compatibilities_to_all={SequenceID('seq0'): CompatibilityToPath(1.0, 1),
                                                  SequenceID('seq1'): CompatibilityToPath(0.9, 1),
                                                  SequenceID('seq2'): CompatibilityToPath(0.95, 1),
                                                  SequenceID('seq3'): CompatibilityToPath(0.6, 1)},
                          consensus_path=[NodeID(0), NodeID(2), NodeID(5), NodeID(6), NodeID(10),
                                          NodeID(12), NodeID(13), NodeID(16)]),
            # no compatibilities to all, no mincomp
            ConsensusNode(parent_node_id=ConsensusNodeID(0),
                          sequences_ids=[SequenceID('seq0'), SequenceID('seq1'), SequenceID('seq2')],
                          consensus_id=ConsensusNodeID(1),
                          consensus_path=[NodeID(0), NodeID(2), NodeID(3), NodeID(6), NodeID(10),
                                          NodeID(11), NodeID(13), NodeID(17)])
        ]

        actual_consensuses_fasta_content = fasta.consensuses_tree_to_fasta(self.pangraph, consensuses_tree)
        expected_consensuses_fasta_content = PangraphTests.get_file_content_as_str(expected_consensuses_fasta_path)
        self.assertEqual(expected_consensuses_fasta_content, actual_consensuses_fasta_content)

if __name__ == '__main__':
    unittest.main()
