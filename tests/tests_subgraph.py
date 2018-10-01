import unittest
from ddt import ddt

from context import Subgraph
from context import Pangraph
from context import Node
from context import nucleotides as n

@ddt
class SubgraphTest(unittest.TestCase):

    def setUp(self):
        nodes = [Node(id=0, base=n.code('A'), in_nodes=[], aligned_to=None),
                 Node(id=1, base=n.code('G'), in_nodes=[], aligned_to=None),
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
        self.expected_pangraph = Pangraph()
        self.expected_pangraph.update_nodes(nodes)
        self.expected_pangraph.set_paths(paths_to_node_ids)

    def test_keep_all_sequences(self):


                ()

        subgraph = subgraph.get_subgraph(graph=graph, pm=pm, paths_to_keep=['testseq0', 'testseq1', 'testseq2', 'testseq3'])


if __name__ == '__main__':
    unittest.main()
