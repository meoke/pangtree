import unittest
from ddt import ddt, data, unpack
from pathlib import Path
from typing import List

from tests.context import Pangraph
from tests.context import Node
from tests.context import pathtools
from tests.context import MultialignmentMetadata
from tests.context import SequenceMetadata
from tests.context import make_base as n
from tests.context import DataType

@ddt
@unittest.skip("po read write must be reimplemented")
class PoWriteReadTest(unittest.TestCase):

    @data(([], ""),
          ([0, 1], "L0L1"))
    @unpack
    def test_get_in_nodes_info(self, in_nodes, expected_info):
        node = Node(node_id=0, base=0, in_nodes=in_nodes, aligned_to=None)
        actual_in_nodes_info = powriter.get_in_nodes_info(node)

        self.assertEqual(actual_in_nodes_info, expected_info)

    @data(([], ""),
          ([0, 1], "S0S1"),
          ([1, 3], "S1S3"))
    @unpack
    def test_get_sources_info(self, sources_ids, expected_info):
        actual_sources_info = powriter.get_sources_info(sources_ids)

        self.assertEqual(actual_sources_info, expected_info)

    @data((None, ""),
          (0, "A0"))
    @unpack
    def test_get_aligned_to_info(self, aligned_to, expected_aligned_to):
        node = Node(node_id=0, base=0, in_nodes=[], aligned_to=aligned_to)
        actual_aligned_to_info = powriter.get_aligned_to_info(node)

        self.assertEqual(actual_aligned_to_info, expected_aligned_to)

