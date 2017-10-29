import unittest

from context import POAGraph
from context import Source, Consensus
from context import Node


class POAGraphTests(unittest.TestCase):
    def test_calculating_compatibility(self):
        p = POAGraph(name='',
                     title='',
                     version='',
                     path='',
                     sources = [Source(currentID=0,
                                        name='source1',
                                        title='source1',
                                        active=True,
                                        nodes_IDs=set([2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15]),
                                        consensusID=1,
                                        weight=0),
                                Source(currentID=1,
                                        name='source2',
                                        title='source2',
                                        active=True,
                                        nodes_IDs=set([0, 1, 3, 4, 6, 7, 8, 12, 13, 14, 16, 17]),
                                        consensusID=0,
                                        weight=0)],
                     consensuses =  [Consensus(currentID=0,
                                                name='CONSENS0',
                                                title='consensus produced by heaviest_bundle, containing 1 seqs',
                                                active=True,
                                                nodes_IDs=set([0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17])),
                                     Consensus(currentID=1,
                                                name='CONSENS1',
                                                title='consensus produced by heaviest_bundle, containing 1 seqs',
                                                active=True,
                                                nodes_IDs=set([2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15]))],
                     nodes = [])
        expected_compatibilities_consensus_0 = [0.8462, 0.9167]
        expected_compatibilities_consensus_1 = [1,0.5833]

        p.calculate_compatibility_to_consensuses()
        self.assertEqual(expected_compatibilities_consensus_0, p.consensuses[0].compatibility_to_sources)
        self.assertEqual(expected_compatibilities_consensus_1, p.consensuses[1].compatibility_to_sources)

    def test_sources_deactivation(self):
        p = POAGraph(name='', title='', version='', path='',
                     sources=[Source(currentID=0,
                                     name='source1',
                                     title='source1',
                                     active=True,
                                     nodes_IDs=set([0,4,5,7,9]),
                                     consensusID=-1,
                                     weight=0),
                              Source(currentID=1,
                                     name='source2',
                                     title='source2',
                                     active=True,
                                     nodes_IDs=set([0,2,4,6]),
                                     consensusID=-1,
                                     weight=0),
                              Source(currentID=2,
                                     name='source3',
                                     title='source3',
                                     active=True,
                                     nodes_IDs=set([1,3,4,6,8]),
                                     consensusID=-1,
                                     weight=0),
                              Source(currentID=-1,
                                     name='source4',
                                     title='source4',
                                     active=False,
                                     nodes_IDs=set([1,2,4,6,7,9]),
                                     consensusID=-1,
                                     weight=0)
                              ],
                     consensuses=[],
                     nodes=[Node(currentID=0, base='A', in_nodes={}, aligned_to=set([1]), sources=set([0,1]), consensuses_count=0),
                            Node(currentID=1, base='C', in_nodes={}, aligned_to=set([0]), sources=set([2]), consensuses_count=0),
                            Node(currentID=2, base='G', in_nodes=set([0]), aligned_to=set([3]), sources=set([1]), consensuses_count=0),
                            Node(currentID=3, base='C', in_nodes=set([1]), aligned_to=set([2]), sources=set([2]), consensuses_count=0),
                            Node(currentID=4, base='G', in_nodes=set([0,2,3]), aligned_to=set([]), sources=set([0,1,2]), consensuses_count=0),
                            Node(currentID=5, base='C', in_nodes=set([4]), aligned_to=set([6]), sources=set([0]), consensuses_count=0),
                            Node(currentID=6, base='A', in_nodes=set([4]), aligned_to=set([5]), sources=set([1,2]), consensuses_count=0),
                            Node(currentID=7, base='T', in_nodes=set([5]), aligned_to=set([8]), sources=set([0]), consensuses_count=0),
                            Node(currentID=8, base='G', in_nodes=set([7]), aligned_to=set([7]), sources=set([2]), consensuses_count=0),
                            Node(currentID=9, base='G', in_nodes=set([6]), aligned_to=set([]), sources=set([0]), consensuses_count=0)])
        p.deactivate_different_then([1,2])
        expected_poagraph = POAGraph(name='', title='', version='', path='',
                     sources=[Source(currentID=-1,
                                     name='source1',
                                     title='source1',
                                     active=False,
                                     nodes_IDs=set([0,4,5,7,9]),
                                     consensusID=-1,
                                     weight=0),
                              Source(currentID=0,
                                     name='source2',
                                     title='source2',
                                     active=True,
                                     nodes_IDs=set([0,2,4,6]),
                                     consensusID=-1,
                                     weight=0),
                              Source(currentID=1,
                                     name='source3',
                                     title='source3',
                                     active=True,
                                     nodes_IDs=set([1,3,4,6,8]),
                                     consensusID=-1,
                                     weight=0),
                              Source(currentID=-1,
                                     name='source4',
                                     title='source4',
                                     active=False,
                                     nodes_IDs=set([1,2,4,6,7,9]),
                                     consensusID=-1,
                                     weight=0)
                              ],
                     consensuses=[],
                     nodes=[Node(currentID=0, base='A', in_nodes=set([]), aligned_to=set([1]), sources=set([1]), consensuses_count=0),
                            Node(currentID=1, base='C', in_nodes=set([]), aligned_to=set([0]), sources=set([2]), consensuses_count=0),
                            Node(currentID=2, base='G', in_nodes=set([0]), aligned_to=set([3]), sources=set([1]), consensuses_count=0),
                            Node(currentID=3, base='C', in_nodes=set([1]), aligned_to=set([2]), sources=set([2]), consensuses_count=0),
                            Node(currentID=4, base='G', in_nodes=set([0,2,3]), aligned_to=set([]), sources=set([1,2]), consensuses_count=0),
                            Node(currentID=-1, base='C', in_nodes=set([4]), aligned_to=set([6]), sources=set([]), consensuses_count=0),
                            Node(currentID=5, base='A', in_nodes=set([4]), aligned_to=set([5]), sources=set([1,2]), consensuses_count=0),
                            Node(currentID=-1, base='T', in_nodes=set([5]), aligned_to=set([8]), sources=set([]), consensuses_count=0),
                            Node(currentID=6, base='G', in_nodes=set([7]), aligned_to=set([7]), sources=set([2]), consensuses_count=0),
                            Node(currentID=-1, base='G', in_nodes=set([6]), aligned_to=set([]), sources=set([]), consensuses_count=0)])
        try:
            self.assertEqual(p, expected_poagraph)
        finally:
            _show_differences(p, expected_poagraph)

    def test_sources_deactivation(self):
        poagraph_with_unassigned_consensuses = POAGraph(name='', title='', version='', path='',
                                     sources=[Source(currentID=-1, #schould become active
                                                     name='source1',
                                                     title='source1',
                                                     active=False,
                                                     nodes_IDs=set([0, 4, 5, 7, 9]),
                                                     consensusID=-1,
                                                     weight=0),
                                              Source(currentID=0, #schould become active
                                                     name='source2',
                                                     title='source2',
                                                     active=True,
                                                     nodes_IDs=set([0, 2, 4, 6]),
                                                     consensusID=-1,
                                                     weight=0),
                                              Source(currentID=1, #should become unactive
                                                     name='source3',
                                                     title='source3',
                                                     active=True,
                                                     nodes_IDs=set([1, 3, 4, 6, 8]),
                                                     consensusID=1,
                                                     weight=0),
                                              Source(currentID=-1, #schould become active
                                                     name='source4',
                                                     title='source4',
                                                     active=False,
                                                     nodes_IDs=set([1, 2, 4, 6, 7, 9]),
                                                     consensusID=-1,
                                                     weight=0)
                                              ],
                                     consensuses=[],
                                     nodes=[Node(currentID=0, base='A', in_nodes=set([]),
                                                                    aligned_to=set([1]), sources=set([1]),
                                                                    consensuses_count=0),
                                           Node(currentID=1, base='C', in_nodes=set([]),
                                                aligned_to=set([0]), sources=set([2]),
                                                consensuses_count=0),
                                           Node(currentID=2, base='G', in_nodes=set([0]),
                                                aligned_to=set([3]), sources=set([1]),
                                                consensuses_count=0),
                                           Node(currentID=3, base='C', in_nodes=set([1]),
                                                aligned_to=set([2]), sources=set([2]),
                                                consensuses_count=0),
                                           Node(currentID=4, base='G', in_nodes=set([0, 2, 3]),
                                                aligned_to=set([]), sources=set([1, 2]),
                                                consensuses_count=0),
                                           Node(currentID=-1, base='C', in_nodes=set([4]),
                                                aligned_to=set([6]), sources=set([]),
                                                consensuses_count=0),
                                           Node(currentID=5, base='A', in_nodes=set([4]),
                                                aligned_to=set([5]), sources=set([1, 2]),
                                                consensuses_count=0),
                                           Node(currentID=-1, base='T', in_nodes=set([5]),
                                                aligned_to=set([8]), sources=set([]),
                                                consensuses_count=0),
                                           Node(currentID=6, base='G', in_nodes=set([7]),
                                                aligned_to=set([7]), sources=set([2]),
                                                consensuses_count=0),
                                           Node(currentID=-1, base='G', in_nodes=set([6]),
                                                aligned_to=set([]), sources=set([]),
                                                consensuses_count=0)])

        poagraph_with_unassigned_consensuses.activate_sources_with_consensus_unassigned()

        expected_poagraph = POAGraph(name='', title='', version='', path='',
                                     sources=[Source(currentID=0,
                                                     name='source1',
                                                     title='source1',
                                                     active=True,
                                                     nodes_IDs=set([0, 4, 5, 7, 9]),
                                                     consensusID=-1,
                                                     weight=0),
                                              Source(currentID=1,
                                                     name='source2',
                                                     title='source2',
                                                     active=True,
                                                     nodes_IDs=set([0, 2, 4, 6]),
                                                     consensusID=-1,
                                                     weight=0),
                                              Source(currentID=-1,
                                                     name='source3',
                                                     title='source3',
                                                     active=False,
                                                     nodes_IDs=set([1, 3, 4, 6, 8]),
                                                     consensusID=1,
                                                     weight=0),
                                              Source(currentID=2,
                                                     name='source4',
                                                     title='source4',
                                                     active=True,
                                                     nodes_IDs=set([1, 2, 4, 6, 7, 9]),
                                                     consensusID=-1,
                                                     weight=0)
                                              ],
                                              consensuses=[],
                                             nodes=[Node(currentID=0, base='A', in_nodes=set([]), aligned_to=set([1]),
                                                         sources=set([0, 1]), consensuses_count=0),
                                                    Node(currentID=1, base='C', in_nodes=set([]), aligned_to=set([0]),
                                                         sources=set([3]), consensuses_count=0),
                                                    Node(currentID=2, base='G', in_nodes=set([0]), aligned_to=set([3]),
                                                         sources=set([1, 3]), consensuses_count=0),
                                                    Node(currentID=-1, base='C', in_nodes=set([1]), aligned_to=set([2]),
                                                         sources=set([]), consensuses_count=0),
                                                    Node(currentID=3, base='G', in_nodes=set([0, 2, 3]), aligned_to=set([]),
                                                         sources=set([0, 1, 3]), consensuses_count=0),
                                                    Node(currentID=4, base='C', in_nodes=set([4]), aligned_to=set([6]),
                                                         sources=set([0]), consensuses_count=0),
                                                    Node(currentID=5, base='A', in_nodes=set([4]), aligned_to=set([5]),
                                                         sources=set([1, 3]), consensuses_count=0),
                                                    Node(currentID=6, base='T', in_nodes=set([5]), aligned_to=set([8]),
                                                         sources=set([0, 3]), consensuses_count=0),
                                                    Node(currentID=-1, base='G', in_nodes=set([7]), aligned_to=set([7]),
                                                         sources=set([]), consensuses_count=0),
                                                    Node(currentID=7, base='G', in_nodes=set([6]), aligned_to=set([]),
                                                         sources=set([0, 3]), consensuses_count=0)])
        try:
            self.assertEqual(poagraph_with_unassigned_consensuses, expected_poagraph)
        finally:
            _show_differences(poagraph_with_unassigned_consensuses, expected_poagraph)

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

if __name__ == '__main__':
    unittest.main()

