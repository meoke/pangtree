import os
import unittest
from pathlib import Path

from ddt import ddt

from tests.context import FromEntrezFastaProvider
from tests.context import pathtools

@ddt
@unittest.skip("slow tests - internet connection required")
class FromEntrezFastaProviderTest_DownloadSequences(unittest.TestCase):

    def test_1_download_AB050936v1_first_eleven_nucleotides(self):
        fasta_provider = FromEntrezFastaProvider(email_address="paulina-ph@wp.pl", use_cache=False)
        sequence_id = "AB050936.1"
        actual_sequence = fasta_provider.get_source(sequence_id, None, None)
        p = pathtools.get_child_path(Path(os.getcwd()), "FastaProviders_Tests/FromEntrezFastaProvider_Tests/AB050936.1.fasta")
        with open(p) as fasta_file_hanlder:
            info_line = fasta_file_hanlder.readline()
            expected_sequence = fasta_file_hanlder.read().upper().replace("\n","")
        self.assertEqual(expected_sequence, actual_sequence)

    def test_2_failed_download(self):
        fasta_provider = FromEntrezFastaProvider(email_address="paulina-ph@wp.pl", use_cache=False)
        sequence_id=""
        with self.assertRaises(Exception) as err:
            actual_sequence = fasta_provider.download_from_ncbi(sequence_id, 0, 11)
            self.assertEqual(str(err), f"Cannot download from Entrez sequence of ID: {sequence_id}")


if __name__ == '__main__':
    unittest.main()
