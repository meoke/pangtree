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
        self.temp_dir = toolkit.create_next_sibling_dir('temp', 'maf_reader_testing')

    def tearDown(self):
        pass
        #toolkit.remove_dir(self.temp_dir)

    @data((None, [range(0,1), range(1,2), range(2, 3), range(3,4)], 4),
          ('all', [range(0,4)], 4),
          ('0:2,3', [range(0,3), range(3,4)], 4))
    @unpack
    def test_get_ranges_for_merge_option(self, input, expected_output, blocks_count):

        actual_output = maf_reader._prepare_merge_ranges(input, blocks_count)

        self.assertEqual(expected_output, actual_output)

    @data(('article',
          'all',
    ['#maf version=1 scoring=roast.v3.3'
    ,'a score=1.0'
    ,'s source1  0 2 + 6 AC-'
    ,'s source2  0 3 + 5 ACT'

    ,'a score=2.0'
    ,'s source1 3 4 + 6 GGTC'
    ,'s source2 4 2 + 5 G-A-'],
    [Source(ID=0, name='source1', title='source1', nodes_count=6, first_node_ID=0, active=True, nodes_IDs=set([0, 1, 3, 4, 5, 7]), consensusID=-1, weight=-1),
    Source(ID=1, name='source2', title='source2', nodes_count=5, first_node_ID=0, active=True, nodes_IDs=set([0, 1, 2, 3, 6]), consensusID=-1, weight=-1)],

    [Node(ID=0, base='A', in_nodes=set(), aligned_to=set(), sources_count = 2, consensuses_count = 0),
    Node(ID=1, base='C', in_nodes={0}, aligned_to=set(), sources_count = 2, consensuses_count = 0),
    Node(ID=2, base='T', in_nodes={1}, aligned_to=(), sources_count = 1, consensuses_count = 0),
    Node(ID=3, base='G', in_nodes={1,2}, aligned_to=(), sources_count = 2, consensuses_count = 0),
    Node(ID=4, base='G', in_nodes={3}, aligned_to=(), sources_count = 1, consensuses_count = 0),
    Node(ID=5, base='T', in_nodes={4}, aligned_to={6}, sources_count = 1, consensuses_count = 0),
    Node(ID=6, base='A', in_nodes={3}, aligned_to={5}, sources_count = 1, consensuses_count = 0),
    Node(ID=7, base='C', in_nodes={6}, aligned_to=(), sources_count = 1, consensuses_count = 0)
    ])
    )
    @unpack
    def test_maf_to_poagraph(self, test_case_name, merge_option, maf_lines, sources, nodes):
        self.maf_path = toolkit.save_text("\n".join(maf_lines), self.temp_dir, 'test.maf')

        poagraphs = maf_reader.parse_to_poagraphs(self.maf_path, merge_option, 'test', toolkit.get_parentdir_name(self.maf_path))
        expected_poagraph = POAGraph(name='test',
                                     title='test0',
                                     version='NOVEMBER',
                                     path='',
                                     sources = sources,
                                     consensuses = [],
                                     nodes = nodes)


        self.assertEqual(expected_poagraph, poagraphs[0])

if __name__ == '__main__':
    unittest.main()
