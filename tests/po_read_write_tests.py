import unittest
import numpy as np
from ddt import ddt, data, unpack

from context import po_reader
from context import po_writer
from context import toolkit
from context import POAGraph
from context import Source, Consensus
from context import Node


@ddt
class PoWriterTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = toolkit.create_next_sibling_dir('files', 'po_writer_testing')
        self.remove_temp_dir=True
        self.po_to_poagraph_testing_data = [
            {'test_name': 'test00 - article',
             'po_file_content': ['VERSION=February',
                                  'NAME=test00',
                                  'TITLE=article',
                                  'LENGTH=18',
                                  'SOURCECOUNT=4',
                                  'SOURCENAME=source1',
                                  'SOURCEINFO=13 2 0 1 source1',
                                  'SOURCENAME=source2',
                                  'SOURCEINFO=12 0 0 0 source2',
                                  'SOURCENAME=CONSENS0',
                                  'SOURCEINFO=15 0 0 0 consensus produced by heaviest_bundle, containing 1 seqs',
                                  'SOURCENAME=CONSENS1',
                                  'SOURCEINFO = 13 0 0 1 consensus produced by heaviest_bundle, containing 1 seqs',
                                  't:S1S2',
                                  'h:L0S1S2',
                                  'p:S0S3',
                                  'k:L1L2S0S1S2S3',
                                  'm:L3S0S1S2S3',
                                  'i:L4S0S2S3A6',
                                  'l:L4S1A5',
                                  'v:L5L6S0S1S2S3',
                                  'r:L7S0S1S2S3',
                                  'p:L8S0S2S3',
                                  'q:L9S0S2S3',
                                  'k:L10S0S2S3',
                                  'n:L8L11S0S1S2S3',
                                  'e:L12S0S1S2S3',
                                  't:L13S0S1S2S3',
                                  'v:L14S0S3A16',
                                  'i:L14S1S2A15',
                                  'm:L16S1S2'],
             'poagraph': POAGraph(  name='test00',
                                    title='article',
                                    version='February',
                                    path=self.temp_dir,
                                    sources=[ Source(ID=0, name='source1', title='source1', weight=0, consensus_ID=1),
                                              Source(ID=1,name='source2', title='source2', weight=0, consensus_ID=0)],
                                    consensuses = [Consensus(ID=0, name='CONSENS0',
                                                            title='consensus produced by heaviest_bundle, containing 1 seqs'),
                                                    Consensus(ID=1, name='CONSENS1',
                                                           title='consensus produced by heaviest_bundle, containing 1 seqs')],
                                    nodes = [Node(ID=0, base='T', in_nodes=np.array([]), aligned_to=None),  # , consensuses_count=1),
                                              Node(ID=1, base='H', in_nodes=np.array([0]), aligned_to=None),  # , consensuses_count=1),
                                              Node(ID=2, base='P', in_nodes=np.array([]), aligned_to=None),  # , consensuses_count=1),
                                              Node(ID=3, base='K', in_nodes=np.array([1, 2]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=4, base='M', in_nodes=np.array([3]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=5, base='I', in_nodes=np.array([4]), aligned_to=6),  # , consensuses_count=2),
                                              Node(ID=6, base='L', in_nodes=np.array([4]), aligned_to=5),
                                              Node(ID=7, base='V', in_nodes=np.array([5, 6]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=8, base='R', in_nodes=np.array([7]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=9, base='P', in_nodes=np.array([8]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=10, base='Q', in_nodes=np.array([9]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=11, base='K', in_nodes=np.array([10]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=12, base='N', in_nodes=np.array([8, 11]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=13, base='E', in_nodes=np.array([12]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=14, base='T', in_nodes=np.array([13]), aligned_to=None),  # , consensuses_count=2),
                                              Node(ID=15, base='V', in_nodes=np.array([14]), aligned_to=16),  # , consensuses_count=1),
                                              Node(ID=16, base='I', in_nodes=np.array([14]), aligned_to=15),  # , consensuses_count=1),
                                              Node(ID=17, base='M', in_nodes=np.array([16]), aligned_to=None)],
                                    ns=np.array([[False, False, True, True, True, True, False, True, True, True, True, True, True, True, True, True,
                                                    False, False],
                                                   [True, True, False, True, True, False, True, True, True, False, False, False, True, True, True,
                                                    False, True, True]]),
                                    nc=np.array([[True, True, False, True, True, True, False, True, True, True, True, True, True, True, True, False,
                                True, True],
                                                    [False, False, True, True, True, True, False, True, True, True, True, True, True, True, True, True,
                                                    False, False]]))},
            {'test_name': 'test01 - 4 nodes aligned',
             'po_file_content': ['VERSION=February',
                                  'NAME=test01',
                                  'TITLE=4 nodes aligned',
                                  'LENGTH=4',
                                  'SOURCECOUNT=5',
                                  'SOURCENAME=source1',
                                  'SOURCEINFO=1 0 0 0 source1',
                                  'SOURCENAME=source2',
                                  'SOURCEINFO=1 0 0 -1 source2',
                                  'SOURCENAME=source3',
                                  'SOURCEINFO=1 0 0 -1 source3',
                                  'SOURCENAME=source4',
                                  'SOURCEINFO=1 0 0 -1 source4',
                                  'SOURCENAME=CONSENS0',
                                  'SOURCEINFO=1 0 0 0 consensus 0',
                                  't:S0S4A1',
                                  'h:S1A2',
                                  'p:S2A3',
                                  'k:S3A0'],
             'poagraph': POAGraph(name='test01',
                                  title='4 nodes aligned',
                                  version='February',
                                  path=self.temp_dir,
                                  sources=[Source(ID=0, name='source1', title='source1', weight=0, consensus_ID=0),
                                          Source(ID=1, name='source2', title='source2', weight=0, consensus_ID=-1),
                                          Source(ID=2, name='source3', title='source3', weight=0, consensus_ID=-1),
                                          Source(ID=3, name='source4', title='source4', weight=0, consensus_ID=-1)],
                                  consensuses=[Consensus(ID=0, name='CONSENS0', title='consensus 0')],
                                  nodes=[Node(ID=0, base='T', in_nodes=np.array([]), aligned_to=1),  # , consensuses_count=0),
                                          Node(ID=1, base='H', in_nodes=np.array([]), aligned_to=2),  # , consensuses_count=0),
                                          Node(ID=2, base='P', in_nodes=np.array([]), aligned_to=3),  # , consensuses_count=0),
                                          Node(ID=3, base='K', in_nodes=np.array([]), aligned_to=0)],
                                  ns=np.array([[True, False, False, False],
                                                           [False, True, False, False],
                                                           [False, False, True, False],
                                                           [False, False, False, True]]),
                                  nc=np.array([True, False, False, False]))},

        ]
        self.poagraph_to_po_testing_data = [
            {'test_name': 'test00 - article',
             'po_file_content': ['VERSION=February',
                                 'NAME=test00',
                                 'TITLE=article',
                                 'LENGTH=18',
                                 'SOURCECOUNT=2',
                                 'SOURCENAME=source1',
                                 'SOURCEINFO=13 2 0 -1 source1',
                                 'SOURCENAME=source2',
                                 'SOURCEINFO=12 0 100 -1 source2',
                                 't:S1',
                                 'h:L0S1',
                                 'p:S0',
                                 'k:L1L2S0S1',
                                 'm:L3S0S1',
                                 'i:L4S0A6',
                                 'l:L4S1A5',
                                 'v:L5L6S0S1',
                                 'r:L7S0S1',
                                 'p:L8S0',
                                 'q:L9S0',
                                 'k:L10S0',
                                 'n:L8L11S0S1',
                                 'e:L12S0S1',
                                 't:L13S0S1',
                                 'v:L14S0A16',
                                 'i:L14S1A15',
                                 'm:L16S1'],
             'poagraph': POAGraph(name='test00',
                                  title='article',
                                  version='February',
                                  path=self.temp_dir,
                                  sources=[Source(ID=0, name='source1', title='source1', weight=0, consensus_ID=1),
                                           Source(ID=1, name='source2', title='source2', weight=0, consensus_ID=0)],
                                  consensuses=[Consensus(ID=0, name='CONSENS0',
                                                         title='consensus produced by heaviest_bundle, containing 1 seqs'),
                                               Consensus(ID=1, name='CONSENS1',
                                                         title='consensus produced by heaviest_bundle, containing 1 seqs')],
                                  nodes=[Node(ID=0, base='T', in_nodes=np.array([]), aligned_to=None),
                                         # , consensuses_count=1),
                                         Node(ID=1, base='H', in_nodes=np.array([0]), aligned_to=None),
                                         # , consensuses_count=1),
                                         Node(ID=2, base='P', in_nodes=np.array([]), aligned_to=None),
                                         # , consensuses_count=1),
                                         Node(ID=3, base='K', in_nodes=np.array([1, 2]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=4, base='M', in_nodes=np.array([3]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=5, base='I', in_nodes=np.array([4]), aligned_to=6),
                                         # , consensuses_count=2),
                                         Node(ID=6, base='L', in_nodes=np.array([4]), aligned_to=5),
                                         Node(ID=7, base='V', in_nodes=np.array([5, 6]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=8, base='R', in_nodes=np.array([7]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=9, base='P', in_nodes=np.array([8]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=10, base='Q', in_nodes=np.array([9]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=11, base='K', in_nodes=np.array([10]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=12, base='N', in_nodes=np.array([8, 11]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=13, base='E', in_nodes=np.array([12]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=14, base='T', in_nodes=np.array([13]), aligned_to=None),
                                         # , consensuses_count=2),
                                         Node(ID=15, base='V', in_nodes=np.array([14]), aligned_to=16),
                                         # , consensuses_count=1),
                                         Node(ID=16, base='I', in_nodes=np.array([14]), aligned_to=15),
                                         # , consensuses_count=1),
                                         Node(ID=17, base='M', in_nodes=np.array([16]), aligned_to=None)],
                                  ns=np.array([[False, False, True, True, True, True, False, True, True, True, True,
                                                True, True, True, True, True,
                                                False, False],
                                               [True, True, False, True, True, False, True, True, True, False, False,
                                                False, True, True, True,
                                                False, True, True]]),
                                  nc=np.array([[True, True, False, True, True, True, False, True, True, True, True,
                                                True, True, True, True, False,
                                                True, True],
                                               [False, False, True, True, True, True, False, True, True, True, True,
                                                True, True, True, True, True,
                                                False, False]])),
             'sources_IDs_to_keep':np.array([0,1])},
            {'test_name': 'test01 - 4 nodes aligned',
             'po_file_content': ['VERSION=February',
                                 'NAME=test01',
                                 'TITLE=4 nodes aligned',
                                 'LENGTH=4',
                                 'SOURCECOUNT=4',
                                 'SOURCENAME=source1',
                                 'SOURCEINFO=1 0 100 -1 source1',
                                 'SOURCENAME=source2',
                                 'SOURCEINFO=1 1 100 -1 source2',
                                 'SOURCENAME=source3',
                                 'SOURCEINFO=1 2 100 -1 source3',
                                 'SOURCENAME=source4',
                                 'SOURCEINFO=1 3 100 -1 source4',
                                 't:S0A1',
                                 'h:S1A2',
                                 'p:S2A3',
                                 'k:S3A0'],
             'poagraph': POAGraph(name='test01',
                                  title='4 nodes aligned',
                                  version='February',
                                  path=self.temp_dir,
                                  sources=[Source(ID=0, name='source1', title='source1', weight=0, consensus_ID=0),
                                           Source(ID=1, name='source2', title='source2', weight=0, consensus_ID=-1),
                                           Source(ID=2, name='source3', title='source3', weight=0, consensus_ID=-1),
                                           Source(ID=3, name='source4', title='source4', weight=0, consensus_ID=-1)],
                                  consensuses=[Consensus(ID=0, name='CONSENS0', title='consensus 0')],
                                  nodes=[Node(ID=0, base='T', in_nodes=np.array([]), aligned_to=1),
                                         # , consensuses_count=0),
                                         Node(ID=1, base='H', in_nodes=np.array([]), aligned_to=2),
                                         # , consensuses_count=0),
                                         Node(ID=2, base='P', in_nodes=np.array([]), aligned_to=3),
                                         # , consensuses_count=0),
                                         Node(ID=3, base='K', in_nodes=np.array([]), aligned_to=0)],
                                  ns=np.array([[True, False, False, False],
                                               [False, True, False, False],
                                               [False, False, True, False],
                                               [False, False, False, True]]),
                                  nc=np.array([True, False, False, False])),
             'sources_IDs_to_keep':np.array([0,1,2,3])},
            {'test_name': 'test00 - deactivate a source',
             'po_file_content': ['VERSION=February',
                                 'NAME=test02',
                                 'TITLE=deactivate a source',
                                 'LENGTH=11',
                                 'SOURCECOUNT=3',
                                 'SOURCENAME=source1',
                                 'SOURCEINFO=6 0 50 -1 source1',
                                 'SOURCENAME=source2',
                                 'SOURCEINFO=3 1 100 -1 source2',
                                 'SOURCENAME=source4',
                                 'SOURCEINFO=3 2 0 -1 source4',
                                 'c:S0A1',
                                 't:S1A2',
                                 'a:S2A0',
                                 't:L0L1S0S1A4',
                                 'g:L2S2A3',
                                 'a:L3S0',
                                 'c:L5S0A7',
                                 'g:L3S1A6',
                                 'c:L6S0',
                                 't:L8S0A10',
                                 'c:L2S2A9'],
             'poagraph': POAGraph(name='test02',
                                  title='deactivate a source',
                                  version='February',
                                  path=self.temp_dir,
                                  sources=[Source(ID=0, name='source1', title='source1', weight=0, consensus_ID=-1),
                                           Source(ID=1, name='source2', title='source2', weight=0, consensus_ID=-1),
                                           Source(ID=2, name='source3', title='source3', weight=0, consensus_ID=-1),
                                           Source(ID=3, name='source4', title='source4', weight=0, consensus_ID=-1)],
                                  nodes=[Node(ID=0, base='C', in_nodes=np.array([]), aligned_to=1),
                                         Node(ID=1, base='T', in_nodes=np.array([]), aligned_to=2),
                                         Node(ID=2, base='G', in_nodes=np.array([]), aligned_to=3),
                                         Node(ID=3, base='A', in_nodes=np.array([]), aligned_to=0),
                                         Node(ID=4, base='T', in_nodes=np.array([0,1]), aligned_to=5),
                                         Node(ID=5, base='G', in_nodes=np.array([2,3]), aligned_to=4),
                                         Node(ID=6, base='A', in_nodes=np.array([4]), aligned_to=7),
                                         Node(ID=7, base='T', in_nodes=np.array([5]), aligned_to=6),
                                         Node(ID=8, base='C', in_nodes=np.array([6]), aligned_to=9),
                                         Node(ID=9, base='G', in_nodes=np.array([4]), aligned_to=10),
                                         Node(ID=10, base='T', in_nodes=np.array([7]), aligned_to=8),
                                         Node(ID=11, base='C', in_nodes=np.array([8]), aligned_to=12),
                                         Node(ID=12, base='G', in_nodes=np.array([10]), aligned_to=11),
                                         Node(ID=13, base='T', in_nodes=np.array([11,12]), aligned_to=14),
                                         Node(ID=14, base='C', in_nodes=np.array([12,3]), aligned_to=13)],
                                  ns=np.array([[True, False, False, False, True, False, True, False, True, False, False, True, False, True, False],
                                               [False, True, False, False, True, False, False, False, False, True, False, False, False, False, False],
                                               [False, False, True, False, False, True, False, True, False, False, True, False, True, False, True],
                                               [False, False, False, True, False, True, False, False, False, False, False, False, False, False, True]]),
                                  nc=np.array([])),
             'sources_IDs_to_keep':np.array([0,1,3])}

        ]

    def tearDown(self):
        if self.remove_temp_dir:
            toolkit.remove_dir(self.temp_dir)

    @data(('c:L0L1L10S0S123A12', [0, 1, 10], [0, 123], [12]),
          ('a:S0', [], [0], []),
          ('t:S12S135A2', [], [12, 135], [2]),
          ('g:L0L1S13', [0,1], [13], []))
    @unpack
    def test_extract_node_parameters(self, line, expected_l, expected_s, expected_a):
        actual_l, actual_s, actual_a = po_reader._extract_node_parameters(line)

        self.assertEqual(expected_l, actual_l)
        self.assertEqual(expected_s, actual_s)
        self.assertEqual(expected_a, actual_a)

    def test_reading_po_to_poagraph(self):
        for test_data in self.po_to_poagraph_testing_data:
            expected_poagraph = test_data["poagraph"]

            self.maf_path = toolkit.save_text("\n".join(test_data["po_file_content"]), self.temp_dir, 'test.po')
            actual_poagraph = po_reader.parse_to_poagraph(str(self.maf_path), output_dir=self.temp_dir)

            try:
                self.assertEqual(expected_poagraph, actual_poagraph)
            except AssertionError as err:
                _show_differences(expected_poagraph, actual_poagraph)
                raise err

    def test_saving_poagraph_to_po(self):
        for test_data in self.poagraph_to_po_testing_data:
            expected_po_file_content = test_data["po_file_content"]

            po_file_path, _ = po_writer.save_as_po(test_data["poagraph"], test_data["sources_IDs_to_keep"])

            with open(po_file_path) as po_file:
                actual_po_file_content = [*map(lambda line : line.strip(), po_file.readlines())]
            try:
                self.assertEqual(expected_po_file_content, actual_po_file_content)
            except AssertionError as err:
                self.remove_temp_dir = False
                raise err

def _show_differences(poagraph1, poagraph2):
    def compare_objects(obj1, obj2, object_name):
        if obj1 != obj2:
            print(object_name, ": ", str(obj1), str(obj2))

    def compare_sequences(seq1, seq2, sequence_name):
        if len(seq1) != len(seq2):
            print(sequence_name, " have different lengths.")
        for i, obj in enumerate(seq1):
            if obj != seq2[i]:
                print(sequence_name, str(i), '\n', str(obj), '\n', str(seq2[i]))

    def compare_ns(ns1, ns2, ns_name):
        if len(ns1) != len(ns2):
            print(ns_name, " have different lengths.")
        for i, obj in enumerate(ns1):
            if not np.array_equal(obj, ns2[i]):
                print(ns_name, str(i), '\n', str(obj), '\n', str(ns2[i]))

    compare_objects(poagraph1.name, poagraph2.name, "name")
    compare_objects(poagraph1.title, poagraph2.title, "title")
    compare_objects(poagraph1.version, poagraph2.version, "version")
    compare_sequences(poagraph1.nodes, poagraph2.nodes, "nodes")
    compare_sequences(poagraph1.sources, poagraph2.sources, "sources")
    compare_sequences(poagraph1.consensuses, poagraph2.consensuses, "consensuses")
    compare_ns(poagraph1.ns, poagraph2.ns, "ns")
    compare_ns(poagraph1.nc, poagraph2.nc, "nc")


if __name__ == '__main__':
    unittest.main()
