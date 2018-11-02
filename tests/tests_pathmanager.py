import unittest
from ddt import ddt, data, unpack

from context import PathManager

@ddt
class PanthmanagerTest(unittest.TestCase):

    def setUp(self):
        paths_to_node_ids = {
            'testseq0': [0, 2, 4, 6, 7, 8, 12, 14, 16],
            'testseq1': [1, 2, 5, 6, 7, 9],
            'testseq2': [3, 4, 6, 7, 10, 12, 14, 17],
            'testseq3': [11, 13, 14, 15]
        }
        self.pathmanager = PathManager()
        self.pathmanager.init_from_dict(18, paths_to_node_ids)

    @data((0, [0]),
          (2, [0, 1]),
          (14, [0, 2, 3]))
    @unpack
    def test_get_sources_ids(self, node_id, expected_sources_ids):
        actual_sources_ids = self.pathmanager.get_sources_ids(node_id)
        self.assertListEqual(actual_sources_ids, expected_sources_ids)