import unittest
from pathlib import Path

from ddt import unpack, data, ddt

from tests.context import missings, graph, multialignment


@ddt
class FromNCBITests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_provider = missings.FromNCBI(use_cache=False)

    @unittest.skip("slow test - internet connection required")
    def test_0_get_10th_symbol_of_AB050936v1(self):
        sequence_id = multialignment.SequenceID("AB050936.1", skip_part_before_dot=False)
        actual_base = self.fasta_provider.get_base(sequence_id, 10)
        path = Path('tests/data/fasta_providers/fasta_files/AB050936.1.fasta')
        expected_base = self.read_sequence(path)[10]
        self.assertEqual(expected_base, actual_base)

    @unittest.skip("slow test - internet connection required")
    def test_1_download_AB050936v1(self):
        fasta_provider = missings.FromNCBI(use_cache=False)
        sequence_id = multialignment.SequenceID("AB050936.1", skip_part_before_dot=False)
        actual_sequence = fasta_provider._download_from_ncbi(sequence_id)
        p = Path('tests/data/fasta_providers/fasta_files/AB050936.1.fasta')
        expected_sequence = self.read_sequence(p)
        self.assertEqual(expected_sequence, actual_sequence)

    @unittest.skip("slow test - internet connection required")
    def test_2_failed_download(self):
        fasta_provider = missings.FromNCBI(use_cache=False)
        sequence_id = ""
        with self.assertRaises(Exception) as err:
            _ = fasta_provider._download_from_ncbi(SequenceID(sequence_id))
            self.assertEqual(str(err), f"Cannot download from Entrez sequence of ID: {sequence_id}")

    @data((multialignment.SequenceID("plain", False), "plain"),
          (multialignment.SequenceID("with.dot", False), "with.dot"),
          (multialignment.SequenceID("with.two.dots", False), "with.two.dots"),
          (multialignment.SequenceID("withv1", False), "with.1"))
    @unpack
    def test_3_guess_entrez_id(self, sequenceID: multialignment.SequenceID, expected_guessed_entrez_id: str):
        fasta_provider = missings.FromNCBI(use_cache=False)
        actual_guessed_entrez_id = fasta_provider._guess_ncbi_sequence_id(sequenceID)

        self.assertEqual(expected_guessed_entrez_id, actual_guessed_entrez_id)

    @staticmethod
    def read_sequence(path: Path):
        with open(path) as fasta_file_hanlder:
            _ = fasta_file_hanlder.readline()
            return fasta_file_hanlder.read().upper().replace("\n", "")


if __name__ == '__main__':
    unittest.main()
