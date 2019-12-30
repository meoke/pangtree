import unittest
from pathlib import Path

from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import missings, msa



class FromFileFastaProviderFastaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_dir = Path(__file__).parent.joinpath("files_zip/").resolve()

    @staticmethod
    def read_sequence(path: Path):
        with open(path) as fasta_file_hanlder:
            _ = fasta_file_hanlder.readline()
            return fasta_file_hanlder.read().upper().replace("\n", "")

    def raise_error_if_unequal(self,
                               sequence_id: msa.SequenceID,
                               expected_sequence: str,
                               fasta_provider: missings.FromFile) -> None:
        for i, expected_symbol in enumerate(expected_sequence):
            expected_base = graph.Base(expected_symbol)
            actual_base = fasta_provider.get_base(sequence_id, i)
            self.assertEqual(expected_base, actual_base)

    def test_1_one_sequence_one_file_in_zip(self):
        fasta_path = self.fasta_dir.joinpath("test_1_one_sequence_one_file_in_zip.zip")

        fasta_provider = missings.FromFile(Path(fasta_path))

        sequence_id = msa.SequenceID("seq1")
        expected_sequence = "ACTGGGTGGGA"

        self.raise_error_if_unequal(sequence_id, expected_sequence, fasta_provider)

    def test_2_three_sequences_in_two_files_in_zip(self):
        fasta_path = self.fasta_dir.joinpath("test_2_three_sequences_in_two_files_in_zip.zip")

        fasta_provider = missings.FromFile(Path(fasta_path))

        sequence_id_1 = msa.SequenceID("seq1")
        self.raise_error_if_unequal(sequence_id_1, "ACTGGGTGGGA", fasta_provider)

        sequence_id_2 = msa.SequenceID("seq2")
        self.raise_error_if_unequal(sequence_id_2, "AA", fasta_provider)

        sequence_id_3 = msa.SequenceID("seq3")
        self.raise_error_if_unequal(sequence_id_3, "GT", fasta_provider)

    def test_3_empty_sequence_name(self):
        fasta_path = self.fasta_dir.joinpath("test_3_empty_sequence_name.zip")

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "No sequences in zipped fastas or incorrect zipped files."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_4_empty_sequence(self):
        fasta_path = self.fasta_dir.joinpath("test_4_empty_sequence.zip")

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "Empty sequence in FASTA. Provide the sequence or remove its header."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_5_empty_fasta(self):
        fasta_path = self.fasta_dir.joinpath("test_5_empty_fasta.zip")

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "No sequences in zipped fastas or incorrect zipped files."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_6_empty_zip(self):
        fasta_path = self.fasta_dir.joinpath("test_6_empty_zip.zip")

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "Incorrect zip fasta source."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_7_no_fasta_in_zip(self):
        fasta_path = self.fasta_dir.joinpath("test_7_no_fasta_in_zip.zip")

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "No sequences in zipped fastas or incorrect zipped files."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)


if __name__ == '__main__':
    unittest.main()
