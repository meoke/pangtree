import unittest

from ddt import ddt, data, unpack
from tests.context import MAX1, MAX2, NODE1, NODE2, NODE3, NODE4
import numpy as np

import unittest
from ddt import ddt
import numpy as np
from tests.context import Pangraph
from tests.context import Node
from tests.context import make_nucleobase as n
from tests.context import SequenceID
from tests.context import ConsensusNode
from tests.context import TreePOAConsensusGenerator
from tests.context import MAX1, NODE4
from tests.context import CompatibilityToPath

@ddt
class CutoffDependingOnPParameter_Tests(unittest.TestCase):

    def setUp(self):
        nodes = [Node(node_id=0, base=n('T'), aligned_to=None),
                 Node(node_id=1, base=n('A'), aligned_to=None),
                 Node(node_id=2, base=n('G'), aligned_to=None),
                 Node(node_id=3, base=n('A'), aligned_to=None),
                 Node(node_id=4, base=n('C'), aligned_to=None),
                 Node(node_id=5, base=n('A'), aligned_to=None),
                 Node(node_id=6, base=n('C'), aligned_to=None),

                 Node(node_id=7, base=n('G'), aligned_to=None),
                 Node(node_id=8, base=n('T'), aligned_to=None),
                 Node(node_id=9, base=n('A'), aligned_to=None)
                 ]

        paths = {
            'testseq0': [[10, 11, 12, 13, 14, 15, 16, 17, 18, 9]],
            'testseq1': [[10, 11, 12, 13, 14, 15, 16, 17, 8, 9]],
            'testseq2': [[10, 11, 12, 13, 14, 15, 16, 7, 8, 9]],
            'testseq3': [[10, 11, 12, 3, 4, 5, 6, 7, 8, 9]],
            'testseq4': [[10, 11, 2, 3, 4, 5, 6, 7, 8, 9]]
        }
        self.pangraph = Pangraph()
        self.pangraph.nodes = nodes
        self.pangraph.paths = paths

    @data(
        (0.5, 0.836660026534076
),
        (1, 0.7),
        (4, 0.6561
),
        )
    @unpack
    def test_max2_p_05_1_4(self, p_value, expected_cutoff):

        consensus_path = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        compatibilities = self.pangraph.get_compatibilities(self.pangraph.get_sequences_ids(), consensus_path, p_value)

        max2_strategy: MAX2 = MAX2()
        actual_max2_cutoff = max2_strategy.find_max_cutoff([*compatibilities.values()]).cutoff.value
        self.assertAlmostEqual(expected_cutoff, actual_max2_cutoff)


if __name__ == '__main__':
    unittest.main()
