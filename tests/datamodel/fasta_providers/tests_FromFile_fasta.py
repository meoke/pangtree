import unittest
from pathlib import Path


from tests.context import missings, graph, multialignment


class FromFileFastaProviderFastaTests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_dir = "tests/datamodel/fasta_providers/files_fasta/"

    @staticmethod
    def read_sequence(path: Path):
        with open(path) as fasta_file_hanlder:
            _ = fasta_file_hanlder.readline()
            return fasta_file_hanlder.read().upper().replace("\n", "")

    def raise_error_if_unequal(self,
                               sequence_id: multialignment.SequenceID,
                               expected_sequence: str,
                               fasta_provider: missings.FromFile) -> None:
        for i, expected_symbol in enumerate(expected_sequence):
            expected_base = graph.Base(expected_symbol)
            actual_base = fasta_provider.get_base(sequence_id, i)
            self.assertEqual(graph.Base(expected_symbol), actual_base)

    def test_1_one_sequence(self):
        fasta_path = self.fasta_dir + "test_1_one_sequence.fasta"
        fasta_provider = missings.FromFile(Path(fasta_path))

        sequence_id = multialignment.SequenceID("seq1")
        expected_sequence = self.read_sequence(fasta_path)

        self.raise_error_if_unequal(sequence_id, expected_sequence, fasta_provider)

    def test_2_three_sequences(self):
        fasta_path = self.fasta_dir + "test_2_three_sequences.fasta"

        fasta_provider = missings.FromFile(Path(fasta_path))

        sequence_id_1 = multialignment.SequenceID("seq1")
        self.raise_error_if_unequal(sequence_id_1, "ACTGGGTGGGA", fasta_provider)

        sequence_id_2 = multialignment.SequenceID("seq2")
        self.raise_error_if_unequal(sequence_id_2, "AA", fasta_provider)

        sequence_id_3 = multialignment.SequenceID("seq3")
        self.raise_error_if_unequal(sequence_id_3, "GT", fasta_provider)

    def test_3_empty_sequence_name(self):
        fasta_path = self.fasta_dir + "test_3_empty_sequence_name.fasta"

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "No _sequences in fasta provided as fasta source or incorrect fasta."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_4_empty_sequence(self):
        fasta_path = self.fasta_dir + "test_4_empty_sequence.fasta"

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "Empty sequence in fasta source file. " \
                           "Provide the sequence or remove its identifier."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_5_empty_fasta(self):
        fasta_path = self.fasta_dir + "test_5_empty_fasta.fasta"

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "No _sequences in fasta provided as fasta source or incorrect fasta."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)

    def test_6_repeated_sequences(self):
        fasta_path = self.fasta_dir + "test_6_repeated_sequences.fasta"

        with self.assertRaises(Exception) as exp:
            _ = missings.FromFile(Path(fasta_path))

        expected_message = "Incorrect fasta provided as fasta source. Sequences ids are not unique."
        actual_message = str(exp.exception)
        self.assertEqual(expected_message, actual_message)


if __name__ == '__main__':
    unittest.main()
