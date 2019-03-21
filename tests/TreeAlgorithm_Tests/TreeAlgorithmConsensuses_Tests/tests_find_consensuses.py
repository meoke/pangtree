import os
import unittest
from pathlib import Path

from ddt import ddt
from tests.context import Node, Pangraph
from tests.context import make_nucleobase as n
from tests.context import TreePOAConsensusGenerator, MAX1, MAX2, NODE1, NODE2, NODE3, NODE4
from tests.context import Compatibility

@ddt
class FindConsensusesTests(unittest.TestCase):

    def setUp(self):
        self.blosum_path = Path(os.path.abspath(__file__)).joinpath('../../../bin/blosum80.mat')

    def test1(self):
        nodes = [Node(node_id=0, base=n('A'), aligned_to=1),
                 Node(node_id=1, base=n('C'), aligned_to=2),
                 Node(node_id=2, base=n('T'), aligned_to=0),
                 Node(node_id=3, base=n('G'), aligned_to=None),
                 Node(node_id=4, base=n('T'), aligned_to=None),
                 Node(node_id=5, base=n('A'), aligned_to=6),
                 Node(node_id=6, base=n('C'), aligned_to=7),
                 Node(node_id=7, base=n('G'), aligned_to=8),
                 Node(node_id=8, base=n('T'), aligned_to=5),
                 Node(node_id=9, base=n('A'), aligned_to=10),
                 Node(node_id=10, base=n('G'), aligned_to=9),

                 Node(node_id=11, base=n('C'), aligned_to=12),
                 Node(node_id=12, base=n('T'), aligned_to=11),
                 Node(node_id=13, base=n('A'), aligned_to=14),
                 Node(node_id=14, base=n('C'), aligned_to=13),
                 Node(node_id=15, base=n('C'), aligned_to=16),
                 Node(node_id=16, base=n('G'), aligned_to=15),
                 Node(node_id=17, base=n('C'), aligned_to=None),
                 Node(node_id=18, base=n('A'), aligned_to=19),
                 Node(node_id=19, base=n('C'), aligned_to=20),
                 Node(node_id=20, base=n('G'), aligned_to=21),
                 Node(node_id=21, base=n('T'), aligned_to=18),

                 Node(node_id=22, base=n('G'), aligned_to=23),
                 Node(node_id=23, base=n('T'), aligned_to=22),
                 Node(node_id=24, base=n('A'), aligned_to=25),
                 Node(node_id=25, base=n('C'), aligned_to=26),
                 Node(node_id=26, base=n('G'), aligned_to=24),
                 Node(node_id=27, base=n('A'), aligned_to=28),
                 Node(node_id=28, base=n('C'), aligned_to=29),
                 Node(node_id=29, base=n('T'), aligned_to=27),
                 Node(node_id=30, base=n('C'), aligned_to=31),
                 Node(node_id=31, base=n('G'), aligned_to=30),

                 Node(node_id=32, base=n('A'), aligned_to=33),
                 Node(node_id=33, base=n('G'), aligned_to=34),
                 Node(node_id=34, base=n('T'), aligned_to=32),
                 Node(node_id=35, base=n('A'), aligned_to=36),
                 Node(node_id=36, base=n('C'), aligned_to=37),
                 Node(node_id=37, base=n('T'), aligned_to=35),
                 Node(node_id=38, base=n('C'), aligned_to=39),
                 Node(node_id=39, base=n('G'), aligned_to=40),
                 Node(node_id=40, base=n('T'), aligned_to=38),
                 Node(node_id=41, base=n('C'), aligned_to=None)
                 ]

        paths = {
            'seq0': [[0, 3, 4, 7, 10, 14, 16, 17, 18, 22, 25, 27, 30, 33, 36, 40]],
            'seq1': [[0, 3, 4, 6, 10, 14, 16, 17, 21, 22, 26, 27, 30, 33, 36, 40]],
            'seq2': [[0, 3, 4, 8, 10, 12, 14, 16, 17, 19, 22, 25, 27, 30, 33, 36, 40]],
            'seq3': [[1, 3, 4, 5, 10, 12, 16, 17, 19, 22, 26, 27, 30, 32, 37, 40]],
            'seq4': [[2, 4, 8, 9, 11, 13, 15, 17, 20, 23, 24, 28, 31, 34, 35, 38, 41]],
            'seq5': [[2, 4, 7, 9, 11, 13, 15, 17, 20, 23, 24, 29, 31, 34, 35, 39, 41]]
        }
        self.pangraph = Pangraph()
        self.pangraph.nodes = nodes
        self.pangraph.paths = paths

        consensuses_generator = TreePOAConsensusGenerator(max_node_strategy=MAX2(),
                                                          node_cutoff_strategy=NODE3(),
                                                          stop=Compatibility(0.99),
                                                          re_consensus=True,
                                                          blosum_path=self.blosum_path
                                                          )
        consensuses_generator.get_consensuses_tree(pangraph=self.pangraph,
                                                   output_dir=Path("TreeAlgorithm_Tests/TreeAlgorithmConsensuses_Tests/output"),
                                                   log_tresholds=False)
        print("KONIEC")


if __name__ == '__main__':
    unittest.main()
