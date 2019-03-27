import os
import unittest
from pathlib import Path
import pandas as pd
import unittest as unittest
from ddt import ddt

from tests.context import pathtools
from tests.context import MultialignmentMetadata

@ddt
class MultialignmentMetadataTest_ReadCSV(unittest.TestCase):

    def test_1_correct(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)

        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],'group': ['A','B','B']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_2_no_seqid(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_2_no_seqid.csv"

        csv_content = pathtools.get_file_content(csv_path)
        with self.assertRaises(Exception) as err:
            _ = MultialignmentMetadata(csv_content)
        self.assertEqual(f"No \'seqid\' column in csv metadata.", str(err.exception))

    def test_3_empty_file(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_3_empty_file.csv"

        csv_content = pathtools.get_file_content(csv_path)
        with self.assertRaises(Exception) as err:
            _ = MultialignmentMetadata(csv_content)
        self.assertEqual(f"Empty csv file.", str(err.exception))

    def test_4_seqid_is_last(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_4_seqid_is_last.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)

        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],'group': ['A','B','B']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_5_double_seqid(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_5_double_seqid.csv"

        csv_content = pathtools.get_file_content(csv_path)
        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],
                                                        'seqid.1': ['s21','s22','s23'],
                                                        'group': ['A','B','B']})
        mm = MultialignmentMetadata(csv_content)

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_6_incorrect_commas_number(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_6_incorrect_commas_number.csv"

        csv_content = pathtools.get_file_content(csv_path)
        with self.assertRaises(Exception) as err:
            _ = MultialignmentMetadata(csv_content)
        self.assertEqual(f"CSV metadata error. Different fields number in line 0 than in header line.", str(err.exception))

    def test_7_not_unique_seqids(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_7_not_unique_seqids.csv"

        csv_content = pathtools.get_file_content(csv_path)
        with self.assertRaises(Exception) as err:
            mm = MultialignmentMetadata(csv_content)
        self.assertEqual(f"Not unique values seqid column in metadata file. Make them unique.", str(err.exception))

    def test_8_correct(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_8_seqids_with_dots.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)

        expected_multialignment_df = pd.DataFrame(index=['s1.1','s2.1','s3.10'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],'group': ['A','B','B']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_9_get_seqids(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        expected_seqids = ['s1', 's2', 's3']

        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)

        actual_seqids = mm.get_all_sequences_ids()

        self.assertEqual(expected_seqids, actual_seqids)

if __name__ == '__main__':
    unittest.main()
