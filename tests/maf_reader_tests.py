import unittest
from ddt import ddt, data, unpack

from context import maf_reader
from context import toolkit
from context import POAGraph
from context import Source
from context import Node

@ddt
class MafReaderTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = toolkit.create_next_sibling_dir('files', 'maf_reader_testing')

    def tearDown(self):
        pass
        toolkit.remove_dir(self.temp_dir)

    @data((None, [range(0,1), range(1,2), range(2, 3), range(3,4)], 4),
          ('all', [range(0,4)], 4),
          ('0:2,3', [range(0,3), range(3,4)], 4))
    @unpack
    def test_get_ranges_for_merge_option(self, input, expected_output, blocks_count):

        actual_output = maf_reader._prepare_merge_ranges(input, blocks_count)

        self.assertEqual(expected_output, actual_output)

    @data(
        ('article',
          'all',
    ['#maf version=1 scoring=roast.v3.3'
    ,'a score=1.0'
    ,'s source1 0 2 + 6 AC-'
    ,'s source2 0 3 + 5 ACT'
    , ''
    ,'a score=2.0'
    ,'s source1 3 4 + 6 GGTC'
    ,'s source2 4 2 + 5 G-A-'],
    [Source(currentID=0, name='source1', title='source1', active=True, nodes_IDs=[0, 1, 3, 4, 5, 7], consensusID=-1, weight=-1),
    Source(currentID=1, name='source2', title='source2', active=True, nodes_IDs=[0, 1, 2, 3, 6], consensusID=-1, weight=-1)],

    [Node(currentID=0, base='A', in_nodes=set(), aligned_to=set(), sources = set([0,1]), consensuses_count = 0),
    Node(currentID=1, base='C', in_nodes={0}, aligned_to=set(), sources = set([0,1]), consensuses_count = 0),
    Node(currentID=2, base='T', in_nodes={1}, aligned_to=set(), sources = set([1]), consensuses_count = 0),
    Node(currentID=3, base='G', in_nodes={1,2}, aligned_to=set(), sources = set([0,1]), consensuses_count = 0),
    Node(currentID=4, base='G', in_nodes={3}, aligned_to=set(), sources = set([0]), consensuses_count = 0),
    Node(currentID=5, base='T', in_nodes={4}, aligned_to={6}, sources = set([0]), consensuses_count = 0),
    Node(currentID=6, base='A', in_nodes={3}, aligned_to={5}, sources = set([1]), consensuses_count = 0),
    Node(currentID=7, base='C', in_nodes={5}, aligned_to=set(), sources = set([0]), consensuses_count = 0)
    ]),
        ('empty',
     'all',
     ['#maf version=1 scoring=roast.v3.3'
         , 'a score=1.0'
         , 's source1 0 0 + 0 ---'
         , 's source2 0 0 + 0 ---'],
     [Source(currentID=0, name='source1', title='source1', active=True, nodes_IDs=[], consensusID=-1, weight=-1),
      Source(currentID=1, name='source2', title='source2', active=True, nodes_IDs=[], consensusID=-1, weight=-1)],

     []),
        ('single letter',
         'all',
         ['#maf version=1 scoring=roast.v3.3'
             , 'a score=1.0'
             , 's source1 1 0 + 1 A'
             , 's source2 1 0 + 1 A'
             , 's source3 1 0 + 1 A'
             , 's source4 1 0 + 1 A'],
         [Source(currentID=0, name='source1', title='source1', active=True, nodes_IDs=[0], consensusID=-1, weight=-1),
          Source(currentID=1, name='source2', title='source2', active=True, nodes_IDs=[0], consensusID=-1, weight=-1),
          Source(currentID=2, name='source3', title='source3', active=True, nodes_IDs=[0], consensusID=-1, weight=-1),
          Source(currentID=3, name='source4', title='source4', active=True, nodes_IDs=[0], consensusID=-1, weight=-1)],

         [Node(currentID=0, base='A', in_nodes={}, aligned_to={}, sources=set([0,1,2,3]), consensuses_count=0)]),
        ("not every block contains all sequences",
         "all",
         ["#maf version=1 scoring=roast.v3.3",
          "a score=1",
          "s test.seq0 0 2 + 3 -TG",
          "s test.seq1 0 1 + 3 C--",
          "",
          "a score=1",
          "s test.seq2 0 1 + 4 T",
          "s test.seq0 2 1 + 3 T",
          "",
          "a score=1",
          "s test.seq1 1 2 + 3 T-A-",
          "s test.seq2 1 3 + 4 GA-C"],
         [Source(currentID=0, name='test.seq0', title='test.seq0', active=True, nodes_IDs=[1, 2, 3], consensusID=-1, weight=-1),
          Source(currentID=1, name='test.seq1', title='test.seq1', active=True, nodes_IDs=[0, 4, 7], consensusID=-1, weight=-1),
          Source(currentID=2, name='test.seq2', title='test.seq2', active=True, nodes_IDs=[3, 5, 6, 8], consensusID=-1, weight=-1)
          ],

    [Node(currentID=0, base='C', in_nodes={}, aligned_to={}, sources=set([1]), consensuses_count=0),
    Node(currentID=1, base='T', in_nodes={}, aligned_to={}, sources=set([0]), consensuses_count=0),
    Node(currentID=2, base='G', in_nodes={1},aligned_to={}, sources=set([0]), consensuses_count=0),
    Node(currentID=3, base='T', in_nodes={2},aligned_to={}, sources=set([0,2]), consensuses_count=0),
    Node(currentID=4, base='T', in_nodes={0},aligned_to={5},sources=set([1]), consensuses_count=0),
    Node(currentID=5, base='G', in_nodes={3},aligned_to={4},sources=set([2]), consensuses_count=0),
    Node(currentID=6, base='A', in_nodes={5},aligned_to={}, sources=set([2]), consensuses_count=0),
    Node(currentID=7, base='A', in_nodes={4},aligned_to={}, sources=set([1]), consensuses_count=0),
    Node(currentID=8, base='C', in_nodes={6},aligned_to={}, sources=set([2]), consensuses_count=0)
          ]
         )
    )
    @unpack
    def test_maf_to_poagraph(self, test_case_name, merge_option, maf_lines, sources, nodes):
        self.maf_path = toolkit.save_text("\n".join(maf_lines), self.temp_dir, 'test.maf')

        poagraphs = maf_reader.parse_to_poagraphs(str(self.maf_path), merge_option, 'test', toolkit.get_parentdir_name(self.maf_path))
        expected_poagraph = POAGraph(name='test',
                                     title='test_0',
                                     version='NOVEMBER',
                                     path='',
                                     sources = sources,
                                     consensuses = [],
                                     nodes = nodes)

        try:
            self.assertEqual(expected_poagraph, poagraphs[0])
        except AssertionError as err:
            _show_differences(expected_poagraph, poagraphs[0])
            raise err

if __name__ == '__main__':
    unittest.main()

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
    compare_objects(poagraph1.name, poagraph2.name, "name")
    compare_objects(poagraph1.title, poagraph2.title, "title")
    compare_objects(poagraph1.version, poagraph2.version, "version")
    compare_sequences(poagraph1.nodes, poagraph2.nodes, "nodes")
    compare_sequences(poagraph1.sources, poagraph2.sources, "sources")
    compare_sequences(poagraph1.consensuses, poagraph2.consensuses, "consensuses")

