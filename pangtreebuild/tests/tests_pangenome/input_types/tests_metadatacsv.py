import unittest
from pathlib import Path

from ddt import unpack, data, ddt

from pangtreebuild.pangenome import graph, builder
from pangtreebuild.pangenome.parameters import missings, msa
from pangtreebuild.tools import pathtools

def sid(x): return msa.SequenceID(x)


def sm(x): return graph.SequenceMetadata(x)


@ddt
class MetadataCSVTests(unittest.TestCase):

    def setUp(self) -> None:
        self.csv_files_dir = Path(__file__).parent.joinpath("csv_files").resolve()
        self.alignment_files_dir = Path(__file__).parent.joinpath("alignment_files").resolve()
        self.fasta_provider = missings.ConstBaseProvider(missings.MissingBase())

    def test_1_correct(self):
        metadata_path = self.csv_files_dir.joinpath("test_1_correct.csv")
        csv_content = pathtools.get_file_content_stringio(metadata_path)

        expected_metadata = {msa.SequenceID('s1'): {'name': 'sequence1', 'group': 'A'},
                             msa.SequenceID('s2'): {'name': 'sequence2', 'group': 'B'},
                             msa.SequenceID('s3'): {'name': 'sequence3', 'group': 'B'}
                             }
        m = msa.MetadataCSV(csv_content, metadata_path)
        actual_metadata = m.metadata

        self.assertEqual(expected_metadata, actual_metadata)

    def test_2_no_seqid(self):
        csv_path = self.csv_files_dir.joinpath("test_2_no_seqid.csv")
        csv_content = pathtools.get_file_content_stringio(csv_path)
        with self.assertRaises(Exception) as err:
            _ = msa.MetadataCSV(csv_content, csv_path)
        self.assertEqual(f"No \'seqid\' column in metadata csv.", str(err.exception))

    def test_3_empty_file(self):
        csv_path = self.csv_files_dir.joinpath("test_3_empty_file.csv")

        csv_content = pathtools.get_file_content_stringio(csv_path)
        with self.assertRaises(Exception) as err:
            _ = msa.MetadataCSV(csv_content, csv_path)
        self.assertEqual(f"Empty csv file.", str(err.exception))

    def test_4_seqid_is_last(self):
        metadata_path = self.csv_files_dir.joinpath("test_4_seqid_is_last.csv")
        csv_content = pathtools.get_file_content_stringio(metadata_path)

        expected_metadata = {msa.SequenceID('s1'): {'name': 'sequence1', 'group': 'A'},
                             msa.SequenceID('s2'): {'name': 'sequence2', 'group': 'B'},
                             msa.SequenceID('s3'): {'name': 'sequence3', 'group': 'B'}
                             }
        m = msa.MetadataCSV(csv_content, metadata_path)
        actual_metadata = m.metadata

        self.assertEqual(expected_metadata, actual_metadata)

    def test_5_double_seqid(self):
        csv_path = self.csv_files_dir.joinpath("test_5_double_seqid.csv")

        csv_content = pathtools.get_file_content_stringio(csv_path)
        with self.assertRaises(Exception) as err:
            _ = msa.MetadataCSV(csv_content, csv_path)
        self.assertEqual("Only one \'seqid\' column in metadata csv is allowed.", str(err.exception))

    def test_6_incorrect_commas_number(self):
        csv_path = self.csv_files_dir.joinpath("test_6_incorrect_commas_number.csv")

        csv_content = pathtools.get_file_content_stringio(csv_path)
        with self.assertRaises(Exception) as err:
            _ = msa.MetadataCSV(csv_content, csv_path)
        self.assertEqual("CSV metadata error. Different number of columns in line 0 than in header line.",
                         str(err.exception))

    def test_7_not_unique_seqids(self):
        csv_path = self.csv_files_dir.joinpath("test_7_not_unique_seqids.csv")

        csv_content = pathtools.get_file_content_stringio(csv_path)
        with self.assertRaises(Exception) as err:
            _ = msa.MetadataCSV(csv_content, csv_path)
        self.assertEqual("Repeated values in seqid column in metadata file. Make them unique.", str(err.exception))

    def test_8_correct_with_dots(self):
        metadata_path = self.csv_files_dir.joinpath("test_8_seqids_with_dots.csv")
        csv_content = pathtools.get_file_content_stringio(metadata_path)

        expected_metadata = {msa.SequenceID('s1'): {'name': 'sequence1', 'group': 'A'},
                             msa.SequenceID('s2'): {'name': 'sequence2', 'group': 'B'},
                             msa.SequenceID('s3'): {'name': 'sequence3', 'group': 'B'}
                             }
        m = msa.MetadataCSV(csv_content, metadata_path)
        actual_metadata = m.metadata

        self.assertEqual(expected_metadata, actual_metadata)

    def test_9_get_seqids(self):
        metadata_path = self.csv_files_dir.joinpath("test_1_correct.csv")
        csv_content = pathtools.get_file_content_stringio(metadata_path)

        expected_seqids = [msa.SequenceID('s1'),
                           msa.SequenceID('s2'),
                           msa.SequenceID('s3')]

        m = msa.MetadataCSV(csv_content, metadata_path)
        actual_seqids = m.get_all_sequences_ids()

        self.assertEqual(expected_seqids, actual_seqids)

    @data(("test_1_the_same_seqid_in_maf_as_in_csv_1",
           "test_10_1.maf", "test_1_correct.csv", "test_10_1.po", {
            sid('s1'): sm({'group': 'A', 'name': 'sequence1'}),
            sid('s2'): sm({'group': 'B', 'name': 'sequence2'}),
            sid('s3'): sm({'group': 'B', 'name': 'sequence3'})}),

          ("test_2_seqid_in_maf_as_in_csv_1_but_with_single_dots",
           "test_10_2.maf", "test_1_correct.csv", "test_10_2.po",
           {sid('s1'): sm({'group': 'A', 'name': 'sequence1'}),
            sid('s2'): sm({'group': 'B', 'name': 'sequence2'}),
            sid('s3'): sm({'group': 'B', 'name': 'sequence3'})}),

          ("test_3_seqid_in_maf_as_in_csv_1_but_with_double_dots",
           "test_10_3.maf", "test_8_seqids_with_dots.csv", "test_10_3.po",
           {sid('s1'): sm({'group': 'A', 'name': 'sequence1'}),
            sid('s2'): sm({'group': 'B', 'name': 'sequence2'}),
            sid('s3'): sm({'group': 'B', 'name': 'sequence3'})}),

          ("test_4_more_seqids_in_maf_than_in_csv",
           "test_10_4.maf", "test_1_correct.csv", "test_10_4.po",
           {sid('s1'): sm({'group': 'A', 'name': 'sequence1'}),
            sid('s2'): sm({'group': 'B', 'name': 'sequence2'}),
            sid('s3'): sm({'group': 'B', 'name': 'sequence3'}),
            sid('s4'): sm({'group': None, 'name': None})}),

          ("test_5_more_seqids_in_csv_than_in_maf",
           "test_10_5.maf", "test_1_correct.csv", "test_10_5.po",
           {sid('s1'): sm({'group': 'A', 'name': 'sequence1'}),
            sid('s2'): sm({'group': 'B', 'name': 'sequence2'}),
            sid('s3'): sm({'group': 'B', 'name': 'sequence3'})}),

          ("test_6_seqids_in_csv_and_maf_all_differ",
           "test_10_6.maf", "test_1_correct.csv", "test_10_6.po",
           {sid('s1'): sm({'group': 'A', 'name': 'sequence1'}),
            sid('s2'): sm({'group': 'B', 'name': 'sequence2'}),
            sid('s3'): sm({'group': 'B', 'name': 'sequence3'}),
            sid('ss1'): sm({'group': None, 'name': None}),
            sid('ss2'): sm({'group': None, 'name': None})})
          )
    @unpack
    def test_10_metadata_feed_to_alignment_from_csv(self,
                                                    test_name,
                                                    maf_name,
                                                    csv_name,
                                                    po_name,
                                                    expected_metadata):
        maf_path = self.alignment_files_dir.joinpath(maf_name)
        csv_path = self.csv_files_dir.joinpath(csv_name)
        po_path = self.alignment_files_dir.joinpath(po_name)

        poagraph, _ = builder.build_from_dagmaf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            self.fasta_provider,
            msa.MetadataCSV(pathtools.get_file_content_stringio(csv_path), csv_path))
        actual_metadata = {seq_id: seq.seqmetadata
                           for seq_id, seq in poagraph.sequences.items()}
        self.assertEqual(expected_metadata, actual_metadata)

        poagraph = builder.build_from_maf(
            msa.Maf(pathtools.get_file_content_stringio(maf_path), maf_path),
            msa.MetadataCSV(pathtools.get_file_content_stringio(csv_path), csv_path))
        actual_metadata = {seq_id: seq.seqmetadata
                           for seq_id, seq in poagraph.sequences.items()}
        self.assertEqual(expected_metadata, actual_metadata)

        poagraph = builder.build_from_po(
            msa.Po(pathtools.get_file_content_stringio(po_path), maf_path),
            msa.MetadataCSV(pathtools.get_file_content_stringio(csv_path), csv_path))
        actual_metadata = {seq_id: seq.seqmetadata
                           for seq_id, seq in poagraph.sequences.items()}
        self.assertEqual(expected_metadata, actual_metadata)


if __name__ == '__main__':
    unittest.main()
