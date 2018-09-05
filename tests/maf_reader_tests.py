import unittest
import numpy as np
from ddt import ddt, data, unpack

from context import maf_reader
from context import toolkit
from context import POAGraph
from context import Source
from context import Node
# from context import Block

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

        actual_output = maf_reader._parse_merge_option_to_ranges(input, blocks_count)

        self.assertEqual(expected_output, actual_output)

    #@unittest.skip("build from maf multialignment")
    @data(
        ('article',
          'all',
    ['#maf version=1 scoring=roast.v3.3'
    ,'a score=1.0'
    ,'s source2 0 3 + 5 ACT'
    ,'s source1 0 2 + 6 AC-'
    , ''
    ,'a score=2.0'
    ,'s source1 3 4 + 6 GGTC'
    ,'s source2 4 2 + 5 G-A-'],
    [Source(ID=0, name='source1', title='source1', weight=0),
    Source(ID=1, name='source2', title='source2', weight=100)],

         [Node(ID=0, base='A', in_nodes=np.array([]), aligned_to=None, consensuses_count=0),
          Node(ID=1, base='C', in_nodes=np.array([0]), aligned_to=None, consensuses_count=0),
          Node(ID=2, base='T', in_nodes=np.array([1]), aligned_to=None, consensuses_count=0),
          Node(ID=3, base='G', in_nodes=np.array([1, 2]), aligned_to=None, consensuses_count=0),
          Node(ID=4, base='G', in_nodes=np.array([3]), aligned_to=None, consensuses_count=0),
          Node(ID=5, base='T', in_nodes=np.array([4]), aligned_to=6, consensuses_count=0),
          Node(ID=6, base='A', in_nodes=np.array([3]), aligned_to=5, consensuses_count=0),
          Node(ID=7, base='C', in_nodes=np.array([5]), aligned_to=None, consensuses_count=0)
          ],
         np.array([[True, True, False, True, True, True, False, True],
                  [True, True, True, True, False, False, True, False]])),
        #TODO CHWILOWO WYRZUCONY TEST
     #    ('empty',
     # 'all',
     # ['#maf version=1 scoring=roast.v3.3'
     #     , 'a score=1.0'
     #     , 's source1 0 0 + 0 ---'
     #     , 's source2 0 0 + 0 ---'],
     # [Source(ID=0, name='source1', title='source1', weight=-1),
     #  Source(ID=1, name='source2', title='source2', weight=-1)],
     #
     # [],
     # np.zeros(shape=(2,0), dtype=np.bool)),
        ('single letter',
         'all',
         ['#maf version=1 scoring=roast.v3.3'
             , 'a score=1.0'
             , 's source1 1 0 + 1 A'
             , 's source2 1 0 + 1 A'
             , 's source3 1 0 + 1 A'
             , 's source4 1 0 + 1 A'],
         [Source(ID=0, name='source1', title='source1',  weight=100),
          Source(ID=1, name='source2', title='source2',  weight=100),
          Source(ID=2, name='source3', title='source3',  weight=100),
          Source(ID=3, name='source4', title='source4',  weight=100)],
         [Node(ID=0, base='A', in_nodes=np.array([]), aligned_to=None, consensuses_count=0)],
         np.array([[True],[True],[True],[True]])),
        ("not every block contains all sequences",
         "all",
         ["#maf version=1 scoring=roast.v3.3",
          "a score=1",
          "s testseq0 0 2 + 3 -TG",
          "s testseq1 0 1 + 3 C--",
          "",
          "a score=1",
          "s testseq2 0 1 + 4 T",
          "s testseq0 2 1 + 3 T",
          "",
          "a score=1",
          "s testseq1 1 2 + 3 T-A-",
          "s testseq2 1 3 + 4 GA-C"],
         [Source(ID=0, name='testseq0', title='testseq0',  weight=100),
          Source(ID=1, name='testseq1', title='testseq1', weight=0),
          Source(ID=2, name='testseq2', title='testseq2',  weight=75)
          ],
         [Node(ID=0, base='C', in_nodes=np.array([]), aligned_to=None, consensuses_count=0),
          Node(ID=1, base='T', in_nodes=np.array([]), aligned_to=None, consensuses_count=0),
          Node(ID=2, base='G', in_nodes=np.array([1]), aligned_to=None, consensuses_count=0),

          Node(ID=3, base='T', in_nodes=np.array([0]), aligned_to=4, consensuses_count=0),
          Node(ID=4, base='G', in_nodes=np.array([8]), aligned_to=3, consensuses_count=0),
          Node(ID=5, base='A', in_nodes=np.array([4]), aligned_to=None, consensuses_count=0),
          Node(ID=6, base='A', in_nodes=np.array([3]), aligned_to=None, consensuses_count=0),
          Node(ID=7, base='C', in_nodes=np.array([5]), aligned_to=None, consensuses_count=0),
          Node(ID=8, base='T', in_nodes=np.array([2]), aligned_to=None, consensuses_count=0),
          ],
         np.array([[False, True, True, False, False, False, False, False, True],
                  [True, False, False, True, False, False, True, False, False],
                  [False, False, False, False, True, True, False, True, True]])
         )
    )
    @unpack
    def test_maf_to_poagraph(self, test_case_name, merge_option, maf_lines, sources, nodes, ns):
        self.maf_path = toolkit.save_text("\n".join(maf_lines), self.temp_dir, 'test.maf')

        poagraphs = maf_reader.parse_to_poagraphs(str(self.maf_path), merge_option, 'test', toolkit.get_parentdir_name(self.maf_path))
        expected_poagraph = POAGraph(name='test',
                                     title='test_0',
                                     version='APRIL',
                                     path='',
                                     sources = sources,
                                     consensuses = [],
                                     nodes = nodes,
                                     ns=ns)

        try:
            self.assertEqual(expected_poagraph, poagraphs[0])
        except AssertionError as err:
            _show_differences(expected_poagraph, poagraphs[0])

            raise err

    @data(
        ('article',
         '',
         ['#maf version=1 scoring=roast.v3.3'
             , 'a score=1.0'
             , 's source1 0 3 + 5 AAA'
             , 's source2 0 3 + 6 BBB'
             , ''
             , 'a score=2.0'
             , 's source1 3 2 + 5 AA-'
             , 's source2 3 2 + 6 B-B'
             , 's source3 0 3 + 5 CCC'
             , ''
             , 'a score=3.0'
             , 's source3 3 2 + 5 CC'
             , 's source2 5 1 + 6 B-'
          ],
         [{0:1, 1:1}, {0:None, 1:2, 2:2}, {1:None, 2:None}]
         )
    )
    @unpack
    @unittest.skip("To implement - Blocks processing")
    def test_maf_to_blocks(self, test_case_name, merge_option, maf_lines, dictionary):
        self.maf_path = toolkit.save_text("\n".join(maf_lines), self.temp_dir, 'test.maf')

        blocks = maf_reader.get_blocks(str(self.maf_path), 'test',
                                                  toolkit.get_parentdir_name(self.maf_path))
        expected_dictionaries = dictionary


        self.assertEqual(expected_dictionaries[0], blocks[0].srcID_to_next_blockID)
        self.assertEqual(expected_dictionaries[1], blocks[1].srcID_to_next_blockID)
        self.assertEqual(expected_dictionaries[2], blocks[2].srcID_to_next_blockID)

if __name__ == '__main__':
    unittest.main()

def _show_differences(poagraph1, poagraph2):
    def compare_objects(obj1, obj2, object_name):
        if obj1 != obj2:
            print(object_name, ": \n", str(obj1), "\n", str(obj2))
    def compare_numpy_arrays(arr1, arr2, name):
        if not np.array_equal(arr1, arr2):
            print(name, ": ", str(arr1), str(arr2))
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
    compare_numpy_arrays(poagraph1.ns, poagraph2.ns, "ns")
    compare_numpy_arrays(poagraph1.nc, poagraph2.ns, "nc")

