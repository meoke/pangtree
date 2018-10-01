import unittest
from ddt import ddt
from pathlib import Path
from textwrap import dedent

from context import Pangraph
from context import Node
from context import nucleotides as n
from context import powriter
from context import pathtools
from context import MultialignmentMetadata
from context import SequenceMetadata

@ddt
class SubgraphTest(unittest.TestCase):

    def setUp(self):
        self.output_dir = pathtools.create_child_dir(Path.cwd(), "po_reader_writer_tests", True)
        genomes_metadata = {0: SequenceMetadata(name="testseq0"),
                            1: SequenceMetadata(name="testseq1"),
                            2: SequenceMetadata(name="testseq2"),
                            3: SequenceMetadata(name="testseq3")}
        self.genomes_info = MultialignmentMetadata(title="test01",
                                                   version="October",
                                                   genomes_metadata=genomes_metadata)
        self.remove_temp_dir = True

    def test_pangraph1(self):
        nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=1),
                 Node(id=1, base=n.code('G'), in_nodes=[], aligned_to=0),
                 Node(id=2, base=n.code('C'), in_nodes=[0, 1], aligned_to=3),
                 Node(id=3, base=n.code('G'), in_nodes=[], aligned_to=2),
                 Node(id=4, base=n.code('A'), in_nodes=[2, 3], aligned_to=5),
                 Node(id=5, base=n.code('T'), in_nodes=[2], aligned_to=4),
                 Node(id=6, base=n.code('G'), in_nodes=[4, 5], aligned_to=None),
                 Node(id=7, base=n.code('G'), in_nodes=[6], aligned_to=None),
                 Node(id=8, base=n.code('A'), in_nodes=[7], aligned_to=9),
                 Node(id=9, base=n.code('C'), in_nodes=[7], aligned_to=10),
                 Node(id=10, base=n.code('G'), in_nodes=[7], aligned_to=11),
                 Node(id=11, base=n.code('T'), in_nodes=[7], aligned_to=8),
                 Node(id=12, base=n.code('A'), in_nodes=[8, 10], aligned_to=13),
                 Node(id=13, base=n.code('C'), in_nodes=[11], aligned_to=12),
                 Node(id=14, base=n.code('T'), in_nodes=[12, 13], aligned_to=None),
                 Node(id=15, base=n.code('A'), in_nodes=[14], aligned_to=16),
                 Node(id=16, base=n.code('C'), in_nodes=[14], aligned_to=17),
                 Node(id=17, base=n.code('G'), in_nodes=[14], aligned_to=15)
                 ]

        paths_to_node_ids = {
            'testseq0': [0, 2, 4, 6, 7, 8, 12, 14, 16],
            'testseq1': [1, 2, 5, 6, 7, 9],
            'testseq2': [3, 4, 6, 7, 10, 12, 14, 17],
            'testseq3': [11, 13, 14, 15]
        }
        pangraph = Pangraph()
        pangraph.update_nodes(nodes)
        pangraph.set_paths(paths_to_node_ids)
        expected_pofile=dedent("""VERSION=September2018
                        NAME=test01_name
                        TITLE=test01_title
                        LENGTH=18
                        SOURCECOUNT=4
                        SOURCENAME=testseq0
                        SOURCEINFO=9 2 0 1 source1
                        SOURCENAME=testseq1
                        SOURCEINFO=6 0 0 0 source2
                        SOURCENAME=testseq2
                        SOURCEINFO=8 x x x d
                        SOURCENAME=testseq3
                        SOURCEINFO=4 0 0 1 coseqs
                        A:S0A1
                        G:S1A0
                        C:L0L1S0S1A3
                        G:S2A2
                        A:L2L3S0S2A5
                        T:L2S1A4
                        G:L4L5S0S1S2
                        G:L7S0S1S2
                        A:L8S0A9
                        C:L7S1A10
                        G:L7S2A11
                        T:L7S3A8
                        A:L8L10S0S2A13
                        C:L11S3A12
                        T:L12L13S0S2S3
                        A:L14S3A16
                        C:L14S0A17
                        G:L14S2A15""")

        #powritertest
        poa_path = pathtools.get_child_file_path(self.output_dir, "pangraph.po")
        powriter.save(pangraph, poa_path, self.genomes_info)
        with open(poa_path) as poa_file:
            actual_pofile = poa_file.read()
        try:
            self.assertEqual(expected_pofile, actual_pofile)
            self.remove_temp_dir = True
        except Exception as ex:
            self.remove_temp_dir = False
            raise ex

        #poreadertest

    def tearDown(self):
        if self.remove_temp_dir:
            pathtools.remove_dir(self.output_dir)


if __name__ == '__main__':
    unittest.main()