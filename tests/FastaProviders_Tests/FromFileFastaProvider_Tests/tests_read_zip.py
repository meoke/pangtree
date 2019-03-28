import unittest
from pathlib import Path

from ddt import ddt

from tests.context import FromFileFastaProvider


@ddt
class FromFileFastaProviderTest_ReadZIP(unittest.TestCase):

    def test_1_one_sequence_one_file_in_zip(self):
        zip_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_1_one_sequence_one_file_in_zip.zip"
        fasta_provider = FromFileFastaProvider(Path(zip_path))

        sequence_id = "seq1"
        actual_sequence = fasta_provider.get_source(sequence_id, None, None)

        expected_sequence = "ACTGGGTGGGA"
        self.assertEqual(expected_sequence, actual_sequence)

    def test_2_three_sequences_in_two_files_in_zip(self):
        zip_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_2_three_sequences_in_two_files_in_zip.zip"
        fasta_provider = FromFileFastaProvider(Path(zip_path))

        sequence_id_1 = "seq1"
        actual_sequence_1 = fasta_provider.get_source(sequence_id_1, None, None)

        sequence_id_2 = "seq2"
        actual_sequence_2 = fasta_provider.get_source(sequence_id_2, None, None)

        sequence_id_3 = "seq3"
        actual_sequence_3 = fasta_provider.get_source(sequence_id_3, None, None)

        expected_sequences = ["ACTGGGTGGGA", "AA", "GT"]
        actual_sequences = [actual_sequence_1, actual_sequence_2, actual_sequence_3]
        self.assertEqual(expected_sequences, actual_sequences)

    def test_3_empty_sequence_name(self):
        fasta_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_3_empty_sequence_name.zip"

        with self.assertRaises(Exception) as exp:
            _ = FromFileFastaProvider(Path(fasta_path))

        expected_message = "No sequences in zip provided as fasta source or incorrect fastas in zip."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_4_empty_sequence(self):
        fasta_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_4_empty_sequence.zip"

        with self.assertRaises(Exception) as exp:
            _ = FromFileFastaProvider(Path(fasta_path))

        expected_message = "Empty sequence in fasta source file. Provide the sequence or remove its identifier."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_5_empty_fasta(self):
        fasta_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_5_empty_fasta.zip"

        with self.assertRaises(Exception) as exp:
            _ = FromFileFastaProvider(Path(fasta_path))

        expected_message = "No sequences in zip provided as fasta source or incorrect fastas in zip."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_6_empty_zip(self):
        fasta_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_6_empty_zip.zip"

        with self.assertRaises(Exception) as exp:
            _ = FromFileFastaProvider(Path(fasta_path))

        expected_message = "Incorrect zip fasta source."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_7_no_fasta_in_zip(self):
        fasta_path = "FastaProviders_Tests/" \
                     "FromFileFastaProvider_Tests/" \
                     "files_zip/" \
                     "test_7_no_fasta_in_zip.zip"

        with self.assertRaises(Exception) as exp:
            _ = FromFileFastaProvider(Path(fasta_path))

        expected_message = "No sequences in zip provided as fasta source or incorrect fastas in zip."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)


if __name__ == '__main__':
    unittest.main()
