import unittest
from typing import List

from ddt import unpack, data, ddt


from pangtreebuild.affinity_tree import builders as at_builders
from pangtreebuild.affinity_tree import parameters as at_params
from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import msa


def sid(x): return msa.SequenceID(x)


def nid(x): return graph.NodeID(x)


def b(x): return graph.Base(x)


@ddt
class AffinityTreeGenerationTests(unittest.TestCase):

    @data((at_params.P(0.5), graph.Compatibility(0.836660026534076)),
          (at_params.P(1), graph.Compatibility(0.7)),
          (at_params.P(4), graph.Compatibility(0.6561)))
    @unpack
    def test_1_p_parameter_influence(self,
                                     p: at_params.P,
                                     expected_cutoff: graph.Compatibility):
        nodes = [graph.Node(node_id=nid(0), base=b('T'), aligned_to=None),
                 graph.Node(node_id=nid(1), base=b('A'), aligned_to=None),
                 graph.Node(node_id=nid(2), base=b('G'), aligned_to=None),
                 graph.Node(node_id=nid(3), base=b('A'), aligned_to=None),
                 graph.Node(node_id=nid(4), base=b('C'), aligned_to=None),
                 graph.Node(node_id=nid(5), base=b('A'), aligned_to=None),
                 graph.Node(node_id=nid(6), base=b('C'), aligned_to=None),
                 graph.Node(node_id=nid(7), base=b('G'), aligned_to=None),
                 graph.Node(node_id=nid(8), base=b('T'), aligned_to=None),
                 graph.Node(node_id=nid(9), base=b('A'), aligned_to=None)
                 ]

        sequences = {
            msa.SequenceID('seq0'):
                graph.Sequence(msa.SequenceID('seq0'),
                               [graph.SeqPath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 18, 9])])],
                               graph.SequenceMetadata({})),
            msa.SequenceID('seq1'):
                graph.Sequence(msa.SequenceID('seq1'),
                               [graph.SeqPath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 8, 9])])],
                               graph.SequenceMetadata({})),
            msa.SequenceID('seq2'):
                graph.Sequence(msa.SequenceID('seq2'),
                               [graph.SeqPath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 7, 8, 9])])],
                               graph.SequenceMetadata({})),
            msa.SequenceID('seq3'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [10, 11, 12, 3, 4, 5, 6, 7, 8, 9])])],
                               graph.SequenceMetadata({})),
            msa.SequenceID('seq4'):
                graph.Sequence(msa.SequenceID('seq3'),
                               [graph.SeqPath([*map(nid, [10, 11, 2, 3, 4, 5, 6, 7, 8, 9])])],
                               graph.SequenceMetadata({}))
        }

        poagraph = graph.Poagraph(nodes, sequences)

        consensus_path = graph.SeqPath([*map(nid, [10, 11, 12, 13, 14, 15, 16, 17, 18, 19])])
        compatibilities = poagraph.get_compatibilities(poagraph.get_sequences_ids(),
                                                       consensus_path,
                                                       p)

        actual_cutoff = at_builders._find_node_cutoff([c for c in compatibilities.values()], []).cutoff
        self.assertAlmostEqual(expected_cutoff.value, actual_cutoff.value)

    @data(
        # single compatibility value
        (0.5, [graph.Compatibility(0.5)]),

        # two compatibilities values
        (0.7, [graph.Compatibility(0.5), graph.Compatibility(0.7)]),
        (1, [graph.Compatibility(1), graph.Compatibility(0.45)]),
        (0.9, [graph.Compatibility(0.9), graph.Compatibility(0.5)]),

        # repeated values
        (0.7, [*map(graph.Compatibility, [0.5, 0.7, 0.7])]),
        (0.9, [*map(graph.Compatibility, [0.9, 0.5, 0.5])]),
        (1, [*map(graph.Compatibility, [0.45, 1, 0.45, 0.45])]),

        # many unique compatibilities values
        (.8, [*map(graph.Compatibility, [.3, .4, .8])]),
        (0.91, [*map(graph.Compatibility, [0.31, 0.32, 0.91, 0.92, 0.93, 0.97])]),
        (0.91, [*map(graph.Compatibility, [0.29, 0.3, 0.33, 0.91, 0.92, 0.93, 0.97])]),
        (1, [*map(graph.Compatibility, [0.81, 0.75, 0.8, 0.81, 1])]),
        (0.9, [*map(graph.Compatibility, [0.5, 0.9, 0.99])]),
        (0.7, [*map(graph.Compatibility, [0.2, 0.85, 0.7, 0.8])]),
        (0.99, [*map(graph.Compatibility, [0.99, 0.9, 0.99])]),
        (0.99, [*map(graph.Compatibility, [0.99])]),

        # repeated distance between values
        (.4, [*map(graph.Compatibility, [.3, .4, .5])]),

        # all the same values
        (.1, [*map(graph.Compatibility, [.1, .1, .1])])
    )
    @unpack
    def test_2_find_cutoff_no_so_far_values(self,
                                            expected_cutoff: float,
                                            compatibilities: List[graph.Compatibility]):
        actual_cutoff = at_builders._find_node_cutoff(compatibilities, []).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)

    def test_3_find_cutoff_no_compatibilities(self):
        with self.assertRaises(ValueError) as err:
            _ = at_builders._find_node_cutoff([], []).cutoff
            self.assertEqual(str(err.exception), """Empty compatibilities list.
                                                    Cannot find cutoff.""")

    @data(
        # guard <= all compatibilities
        (0.2, [0.2, 0.7, 0.8, 0.85], [0.1, 0.01, 0]),
        (0.7, [0.7, 0.85, 0.7, 0.8], [0.1, 0.01, 0]),
        (0.8, [0.7, 0.7, 0.85, 0.8], [0.85, 0.91, 1.0]),

        # guard > all compatibilities
        (0.6, [0.3, 0.6, 0.61, 0.61], [0.99]),  # big distance to guard
        (0.9, [0.2, 0.97, 0.98, 0.9], [0.99]),  # small distance to guard

        # guard between compatibilities
        (0.5, [0.2, 0.57, 0.58, 0.5], [0.55]),  # take smaller than guard
        (0.58, [0.2, 0.27, 0.58, 0.2], [0.55]),  # take greater than guard
        (0.55, [0.2, 0.58, 0.27, 0.55], [0.55])  # take equal to guard
    )
    @unpack
    def test_4_find_cutoff_with_so_far_values(self,
                                              expected_cutoff,
                                              compatibilities,
                                              so_far_cutoffs):
        compatibilities = [graph.Compatibility(c) for c in compatibilities]
        so_far_cutoffs = [graph.Compatibility(c) for c in so_far_cutoffs]
        actual_cutoff = at_builders._find_node_cutoff(compatibilities, so_far_cutoffs).cutoff
        self.assertEqual(expected_cutoff, actual_cutoff.value)


if __name__ == '__main__':
    unittest.main()
