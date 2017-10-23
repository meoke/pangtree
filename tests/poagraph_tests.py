import unittest

from context import POAGraph
from context import Source, Consensus


class POAGraphTests(unittest.TestCase):
    def test_calculating_compatibility(self):
        p = POAGraph(name='',
                     title='',
                     version='',
                     path='',
                     sources = [Source(ID=0,
                                        name='source1',
                                        title='source1',
                                        active=True,
                                        nodes_IDs=set([2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15]),
                                        consensusID=1,
                                        weight=0),
                                Source(ID=1,
                                        name='source2',
                                        title='source2',
                                        active=True,
                                        nodes_IDs=set([0, 1, 3, 4, 6, 7, 8, 12, 13, 14, 16, 17]),
                                        consensusID=0,
                                        weight=0)],
                     consensuses =  [Consensus(ID=0,
                                                name='CONSENS0',
                                                title='consensus produced by heaviest_bundle, containing 1 seqs',
                                                active=True,
                                                nodes_IDs=set([0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17])),
                                     Consensus(ID=1,
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


if __name__ == '__main__':
    unittest.main()
