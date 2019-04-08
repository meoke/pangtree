import unittest
from pathlib import Path

from tests.context import FromNCBI, EmailAddress
from tests.context import pSeq


@unittest.skip("slow test - internet connection required")
class FromNCBITests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_provider = FromNCBI(EmailAddress('a@gmail.com'), use_cache=False)

    def test_0_get_10th_symbol_of_AB050936v1(self):
        sequence_id = pSeq.SequenceID("AB050936.1", skip_part_before_dot=False)
        actual_base = self.fasta_provider.get_base(sequence_id, 10)
        path = Path('tests/data/fasta_providers/fasta_files/AB050936.1.fasta')
        expected_base = self.read_sequence(path)[10]
        self.assertEqual(expected_base, actual_base)

    def test_1_download_AB050936v1(self):
        fasta_provider = FromNCBI(EmailAddress('a@gmail.com'), use_cache=False)
        sequence_id = pSeq.SequenceID("AB050936.1", skip_part_before_dot=False)
        actual_sequence = fasta_provider._download_from_ncbi(sequence_id)
        p = Path('tests/data/fasta_providers/fasta_files/AB050936.1.fasta')
        expected_sequence = self.read_sequence(p)
        self.assertEqual(expected_sequence, actual_sequence)

    def test_2_failed_download(self):
        fasta_provider = FromNCBI(EmailAddress('a@gmail.com'), use_cache=False)
        sequence_id = ""
        with self.assertRaises(Exception) as err:
            _ = fasta_provider._download_from_ncbi(sequence_id)
            self.assertEqual(str(err), f"Cannot download from Entrez sequence of ID: {sequence_id}")

    def read_sequence(self, path: Path):
        with open(path) as fasta_file_hanlder:
            _ = fasta_file_hanlder.readline()
            return fasta_file_hanlder.read().upper().replace("\n", "")


if __name__ == '__main__':
    unittest.main()
