import unittest
import numpy as np
from ddt import ddt, data, unpack

from context import po_reader
from context import toolkit
from context import POAGraph
from context import Source, Consensus
from context import Node

@ddt
class MafReaderTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = toolkit.create_next_sibling_dir('files', 'po_reader_testing')

    def tearDown(self):
        toolkit.remove_dir(self.temp_dir)

    @data(
        ('article',
         ['VERSION=NOVEMBER',
            'NAME=test',
            'TITLE=test_0',
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
            'T:S1S2',
            'H:L0S1S2',
            'P:S0S3',
            'K:L1L2S0S1S2S3',
            'M:L3S0S1S2S3',
            'I:L4S0S2S3A6',
            'L:L4S1A5',
            'V:L5L6S0S1S2S3',
            'R:L7S0S1S2S3',
            'P:L8S0S2S3',
            'Q:L9S0S2S3',
            'K:L10S0S2S3',
            'N:L8L11S0S1S2S3',
            'E:L12S0S1S2S3',
            'T:L13S0S1S2S3',
            'V:L14S0S3A16',
            'I:L14S1S2A15',
            'M:L16S1S2'],
    [Source(ID=0,
            name='source1',
            title='source1',
            weight=0,
            consensus_ID=1),
    Source(ID=1,
           name='source2',
           title='source2',
           weight=0,
           consensus_ID=0)],
    [Consensus(ID=0, 
               name='CONSENS0', 
               title='consensus produced by heaviest_bundle, containing 1 seqs'),
    Consensus(ID=1, 
              name='CONSENS1', 
              title='consensus produced by heaviest_bundle, containing 1 seqs')],
         
    [Node(ID=0,base='T', in_nodes=np.array([]), aligned_to=None),#, consensuses_count=1),
    Node(ID=1, base='H', in_nodes=np.array([0]), aligned_to=None),#, consensuses_count=1),
    Node(ID=2, base='P',in_nodes=np.array([]), aligned_to=None),#, consensuses_count=1),
    Node(ID=3, base='K',in_nodes=np.array([1, 2]), aligned_to=None),#, consensuses_count=2),
    Node(ID=4, base='M',in_nodes=np.array([3]), aligned_to=None),#, consensuses_count=2),
    Node(ID=5, base='I',in_nodes=np.array([4]), aligned_to=6),#, consensuses_count=2),
    Node(ID=6, base='L',in_nodes=np.array([4]), aligned_to=5),
    Node(ID=7, base='V',in_nodes=np.array([5,6]), aligned_to=None),#, consensuses_count=2),
    Node(ID=8, base='R',in_nodes=np.array([7]), aligned_to=None),#, consensuses_count=2),
    Node(ID=9, base='P',in_nodes=np.array([8]), aligned_to=None),#, consensuses_count=2),
    Node(ID=10,base='Q', in_nodes=np.array([9]), aligned_to=None),#, consensuses_count=2),
    Node(ID=11,base='K', in_nodes=np.array([10]), aligned_to=None),#, consensuses_count=2),
    Node(ID=12,base='N', in_nodes=np.array([8,11]), aligned_to=None),#, consensuses_count=2),
    Node(ID=13,base='E', in_nodes=np.array([12]), aligned_to=None),#, consensuses_count=2),
    Node(ID=14,base='T', in_nodes=np.array([13]), aligned_to=None),#, consensuses_count=2),
    Node(ID=15,base='V', in_nodes=np.array([14]), aligned_to=16),#, consensuses_count=1),
    Node(ID=16,base='I', in_nodes=np.array([14]), aligned_to=15),#, consensuses_count=1),
    Node(ID=17,base='M', in_nodes=np.array([16]), aligned_to=None)],
         np.array([[False, False, True, True, True, True, False, True, True, True, True, True, True, True, True, True, False, False],
                   [True, True, False, True, True, False, True, True, True, False, False, False, True, True, True, False, True, True]]),
        np.array([[True, True, False, True, True, True, False, True, True, True, True, True, True, True, True, False, True, True],
                   [False, False, True, True, True, True, False, True, True, True, True, True, True, True, True, True, False, False]
                   ])), #, consensuses_count=1)]),
        ('4 nodes aligned',
         ['VERSION=NOVEMBER',
          'NAME=test',
          'TITLE=test_0',
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
          'T:S0S4A1',
          'H:S1A2',
          'P:S2A3',
          'K:S3A0'],
         [Source(ID=0,
                 name='source1',
                 title='source1',
                 weight=0,
                 consensus_ID=0),
          Source(ID=1,
                 name='source2',
                 title='source2',
                 weight=0,
                 consensus_ID=-1),
         Source(ID=2,
                name='source3',
                title='source3',
                weight=0,
                consensus_ID=-1),
         Source(ID=3,
                name='source4',
                title='source4',
                weight=0,
                consensus_ID=-1)],
         [Consensus(ID=0,
                    name='CONSENS0',
                    title='consensus 0')],
         [Node(ID=0, base='T', in_nodes=np.array([]), aligned_to=1),#, consensuses_count=0),
          Node(ID=1, base='H', in_nodes=np.array([]), aligned_to=2),#, consensuses_count=0),
          Node(ID=2, base='P', in_nodes=np.array([]), aligned_to=3),#, consensuses_count=0),
          Node(ID=3, base='K', in_nodes=np.array([]), aligned_to=0)],
         np.array([[True, False, False, False],
                   [False, True, False, False],
                   [False, False, True, False],
                   [False, False, False, True]]),
         np.array([True, False, False, False]))
    )
    @unpack
    def test_maf_to_poagraph(self, test_case_name, po_lines, expected_sources, expected_consensuses, expected_nodes, ns, nc):
        self.maf_path = toolkit.save_text("\n".join(po_lines), self.temp_dir, 'test.po')

        poagraph = po_reader.parse_to_poagraph(str(self.maf_path), output_dir=self.temp_dir)
        expected_poagraph = POAGraph(name='test',
                                     title='test_0',
                                     version='NOVEMBER',
                                     path='',
                                     sources=expected_sources,
                                     consensuses=expected_consensuses,
                                     nodes=expected_nodes,
                                     ns=ns,
                                     nc=nc)

        try:
            self.assertEqual(expected_poagraph, poagraph)
        except AssertionError as err:
            _show_differences(expected_poagraph, poagraph)
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
