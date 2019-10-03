import unittest
from pathlib import Path


from tests.context import fasta_providers, pSeq, pNode


class FromFileFastaProviderFastaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_dir = "tests/datamodel/fasta_providers/files_zip/"

    @staticmethod
    def read_sequence(path: Path):
        with open(path) as fasta_file_hanlder:
            _ = fasta_file_hanlder.readline()
            return fasta_file_hanlder.read().upper().replace("\n", "")

    def raise_error_if_unequal(self,
                               sequence_id: pSeq.SequenceID,
                               expected_sequence: str,
                               fasta_provider: fasta_providers.FromFile) -> None:
        for i, expected_symbol in enumerate(expected_sequence):
            expected_base = pNode.Base(expected_symbol)
            actual_base = fasta_provider.get_base(sequence_id, i)
            self.assertEqual(expected_base, actual_base)

    def test_1_one_sequence_one_file_in_zip(self):
        fasta_path = self.fasta_dir + "test_1_one_sequence_one_file_in_zip.zip"
        fasta_provider = fasta_providers.FromFile(Path(fasta_path))

        sequence_id = pSeq.SequenceID("seq1")
        expected_sequence = "ACTGGGTGGGA"

        self.raise_error_if_unequal(sequence_id, expected_sequence, fasta_provider)

    def test_2_three_sequences_in_two_files_in_zip(self):
        fasta_path = self.fasta_dir + "test_2_three_sequences_in_two_files_in_zip.zip"

        fasta_provider = fasta_providers.FromFile(Path(fasta_path))

        sequence_id_1 = pSeq.SequenceID("seq1")
        self.raise_error_if_unequal(sequence_id_1, "ACTGGGTGGGA", fasta_provider)

        sequence_id_2 = pSeq.SequenceID("seq2")
        self.raise_error_if_unequal(sequence_id_2, "AA", fasta_provider)

        sequence_id_3 = pSeq.SequenceID("seq3")
        self.raise_error_if_unequal(sequence_id_3, "GT", fasta_provider)

    def test_3_empty_sequence_name(self):
        fasta_path = self.fasta_dir + "test_3_empty_sequence_name.zip"

        with self.assertRaises(Exception) as exp:
            _ = fasta_providers.FromFile(Path(fasta_path))

        expected_message = "No _sequences in zip provided as fasta source or incorrect fastas in zip."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_4_empty_sequence(self):
        fasta_path = self.fasta_dir + "test_4_empty_sequence.zip"

        with self.assertRaises(Exception) as exp:
            _ = fasta_providers.FromFile(Path(fasta_path))

        expected_message = "Empty sequence in fasta source file. " \
                           "Provide the sequence or remove its identifier."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_5_empty_fasta(self):
        fasta_path = self.fasta_dir + "test_5_empty_fasta.zip"

        with self.assertRaises(Exception) as exp:
            _ = fasta_providers.FromFile(Path(fasta_path))

        expected_message = "No _sequences in zip provided as fasta source or incorrect fastas in zip."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_6_empty_zip(self):
        fasta_path = self.fasta_dir + "test_6_empty_zip.zip"

        with self.assertRaises(Exception) as exp:
            _ = fasta_providers.FromFile(Path(fasta_path))

        expected_message = "Incorrect zip fasta source."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_7_no_fasta_in_zip(self):
        fasta_path = self.fasta_dir + "test_7_no_fasta_in_zip.zip"

        with self.assertRaises(Exception) as exp:
            _ = fasta_providers.FromFile(Path(fasta_path))

        expected_message = "No _sequences in zip provided as fasta source or incorrect fastas in zip."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)


if __name__ == '__main__':
    unittest.main()
