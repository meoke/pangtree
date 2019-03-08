import unittest
from pathlib import Path

from ddt import ddt
from tests.context import Node, Pangraph
from tests.context import nucleotides as n
from tests.context import TreePOAConsensusGenerator, MAX1, MAX2, NODE1, NODE2, NODE3, NODE4
from tests.context import Compatibility

@ddt
class FindConsensusesTests(unittest.TestCase):

    def test1(self):
        nodes = [Node(id=0, base=n.code('A'), aligned_to=1),
                 Node(id=1, base=n.code('C'), aligned_to=2),
                 Node(id=2, base=n.code('T'), aligned_to=0),
                 Node(id=3, base=n.code('G'), aligned_to=None),
                 Node(id=4, base=n.code('T'), aligned_to=None),
                 Node(id=5, base=n.code('A'), aligned_to=6),
                 Node(id=6, base=n.code('C'), aligned_to=7),
                 Node(id=7, base=n.code('G'), aligned_to=8),
                 Node(id=8, base=n.code('T'), aligned_to=5),
                 Node(id=9, base=n.code('A'), aligned_to=10),
                 Node(id=10, base=n.code('G'),aligned_to=9),

                 Node(id=11, base=n.code('C'), aligned_to=12),
                 Node(id=12, base=n.code('T'), aligned_to=11),
                 Node(id=13, base=n.code('A'), aligned_to=14),
                 Node(id=14, base=n.code('C'), aligned_to=13),
                 Node(id=15, base=n.code('C'), aligned_to=16),
                 Node(id=16, base=n.code('G'), aligned_to=15),
                 Node(id=17, base=n.code('C'), aligned_to=None),
                 Node(id=18, base=n.code('A'), aligned_to=19),
                 Node(id=19, base=n.code('C'), aligned_to=20),
                 Node(id=20, base=n.code('G'), aligned_to=21),
                 Node(id=21, base=n.code('T'), aligned_to=18),

                 Node(id=22, base=n.code('G'), aligned_to=23),
                 Node(id=23, base=n.code('T'), aligned_to=22),
                 Node(id=24, base=n.code('A'), aligned_to=25),
                 Node(id=25, base=n.code('C'), aligned_to=26),
                 Node(id=26, base=n.code('G'), aligned_to=24),
                 Node(id=27, base=n.code('A'), aligned_to=28),
                 Node(id=28, base=n.code('C'), aligned_to=29),
                 Node(id=29, base=n.code('T'), aligned_to=27),
                 Node(id=30, base=n.code('C'), aligned_to=31),
                 Node(id=31, base=n.code('G'), aligned_to=30),

                 Node(id=32, base=n.code('A'), aligned_to=33),
                 Node(id=33, base=n.code('G'), aligned_to=34),
                 Node(id=34, base=n.code('T'), aligned_to=32),
                 Node(id=35, base=n.code('A'), aligned_to=36),
                 Node(id=36, base=n.code('C'), aligned_to=37),
                 Node(id=37, base=n.code('T'), aligned_to=35),
                 Node(id=38, base=n.code('C'), aligned_to=39),
                 Node(id=39, base=n.code('G'), aligned_to=40),
                 Node(id=40, base=n.code('T'), aligned_to=38),
                 Node(id=41, base=n.code('C'), aligned_to=None)
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
                                                          re_consensus=True
                                                          )
        consensuses_generator.get_consensuses_tree(pangraph=self.pangraph,
                                                   output_dir=Path("TreeAlgorithm_Tests/TreeAlgorithmConsensuses_Tests/output"))
        print("KONIEC")


if __name__ == '__main__':
    unittest.main()
