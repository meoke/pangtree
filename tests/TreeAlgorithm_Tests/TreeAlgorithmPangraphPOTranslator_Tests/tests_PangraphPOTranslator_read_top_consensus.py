import unittest
from ddt import ddt
import numpy as np
from tests.context import PangraphPO_Translator
from tests.context import Pangraph
from tests.context import Node
from tests.context import nucleotides as n
from tests.context import NodePO, SequencePO

@ddt
class PangraphPOTranslator_read_top_consensus_Test(unittest.TestCase):

    def setUp(self):
        nodes = [Node(id=0, base=n.code('T'), aligned_to=None,),
                 Node(id=1, base=n.code('A'), aligned_to=2),
                 Node(id=2, base=n.code('G'), aligned_to=1),
                 Node(id=3, base=n.code('A'), aligned_to=4),
                 Node(id=4, base=n.code('C'), aligned_to=3),
                 Node(id=5, base=n.code('A'), aligned_to=6),
                 Node(id=6, base=n.code('C'), aligned_to=7),
                 Node(id=7, base=n.code('G'), aligned_to=8),
                 Node(id=8, base=n.code('T'), aligned_to=5),
                 Node(id=9, base=n.code('A'), aligned_to=None),
                 Node(id=10, base=n.code('C'), aligned_to=11),
                 Node(id=11, base=n.code('T'), aligned_to=10),
                 Node(id=12, base=n.code('G'), aligned_to=None),
                 Node(id=13, base=n.code('A'), aligned_to=14),
                 Node(id=14, base=n.code('C'), aligned_to=13)
                 ]

        paths = {
            'testseq0': [[0, 1, 3, 5, 9, 10, 13]],
            'testseq1': [[1, 3, 6, 9, 11]],
            'testseq2': [[2, 4, 7, 9, 11, 12]],
            'testseq3': [[2, 4, 8, 9, 11, 12, 14]]
        }
        self.pangraph = Pangraph()
        self.pangraph.nodes = nodes
        self.pangraph.paths = paths

    def test_seq1_only_in_input(self):
        translator = PangraphPO_Translator(self.pangraph, ['testseq1'])
        _ = translator.get_input_po_content()

        poa_output =["VERSION=pangenome\n",
                     "NAME=pangenome\n",
                     "TITLE=pangenome\n",
                     "LENGTH=5\n",
                     "SOURCECOUNT=2\n",
                     "SOURCENAME=testseq1\n",
                     "SOURCEINFO=5 0 100 0 testseq1\n",
                     "SOURCENAME=CONSENS0\n",
                     "SOURCEINFO=5 0 100 0 CONSENS0\n",
                     "a:S0S1\n",
                     "a:L0S0S1\n",
                     "c:L1S0S1\n",
                     "a:L2S0S1\n",
                     "t:L2S1S1"]
        actual_consensus_path = translator.read_top_consensus(poa_output)
        expected_consensus_path = [1, 3, 6, 9, 11]
        self.assertEqual(expected_consensus_path, actual_consensus_path)
