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

    def test_1_the_same_seqid_in_maf_as_in_csv_1(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        mm.feed_with_maf_data(['s1', 's2', 's3'])
        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],
                                                        'group': ['A','B','B'],
                                                        'mafname': ['s1', 's2', 's3']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_2_seqid_in_maf_as_in_csv_1_but_with_single_dots(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        mm.feed_with_maf_data(['db.s1', 'foo.s2', 'foo2.s3'])
        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],
                                                        'group': ['A','B','B'],
                                                        'mafname': ['db.s1', 'foo.s2', 'foo2.s3']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_3_seqid_in_maf_as_in_csv_1_but_with_double_dots(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_8_seqids_with_dots.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        mm.feed_with_maf_data(['db.s1.1', 'foo.s2.1', 'foo2.s3.10'])
        expected_multialignment_df = pd.DataFrame(index=['s1.1','s2.1','s3.10'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],
                                                        'group': ['A','B','B'],
                                                        'mafname': ['db.s1.1', 'foo.s2.1', 'foo2.s3.10']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))


    def test_4_more_seqids_in_maf_than_in_csv(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        mm.feed_with_maf_data(['s1'])
        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3'],
                                                  data={'name': ['sequence1','sequence2','sequence3'],
                                                        'group': ['A','B','B'],
                                                        'mafname': ['s1', None, None]})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))

    def test_5_more_seqids_in_csv_than_in_maf(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        mm.feed_with_maf_data(['s1', 's2', 's3', 's4'])
        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3', 's4'],
                                                  data={'name': ['sequence1','sequence2','sequence3',None],
                                                        'group': ['A','B','B',None],
                                                        'mafname': ['s1', 's2', 's3','s4']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))


    def test_6_seqids_in_csv_and_maf_all_differ(self):
        csv_path = "MultialignmentMetadata_Tests/" \
                   "files_read_csv/" \
                   "test_1_correct.csv"
        csv_content = pathtools.get_file_content(csv_path)
        mm = MultialignmentMetadata(csv_content)
        mm.feed_with_maf_data(['a1', 'a2', 'a3', 'a4'])
        expected_multialignment_df = pd.DataFrame(index=['s1','s2','s3', 'a1','a2','a3','a4'],
                                                  data={'name': ['sequence1','sequence2','sequence3',None,None,None,None],
                                                        'group': ['A','B','B',None,None,None,None],
                                                        'mafname': [None, None, None,'a1', 'a2', 'a3', 'a4']})

        self.assertTrue(expected_multialignment_df.equals(mm.metadata_df))



if __name__ == '__main__':
    unittest.main()
