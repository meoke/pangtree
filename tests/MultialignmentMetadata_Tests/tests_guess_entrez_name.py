import os
import unittest
from pathlib import Path
import pandas as pd
import unittest as unittest
from ddt import ddt, unpack, data

from tests.context import pathtools
from tests.context import MultialignmentMetadata

@ddt
class MultialignmentMetadataTest_ReadCSV(unittest.TestCase):

    @data(("plain", "plain"),
          ("with.dot", "with.dot"),
          ("with.two.dots","with.two.dots"),
          ("withv1","with.1"))
    @unpack
    def test_1(self, sequenceID, expected_guessed_entrez_id):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        actual_guessed_entrez_id = mm.get_entrez_name(sequenceID)

        self.assertEqual(expected_guessed_entrez_id, actual_guessed_entrez_id)


if __name__ == '__main__':
    unittest.main()
