import os
import shutil
import unittest
from pathlib import Path

from ddt import ddt

from tests.context import FromEntrezFastaProvider
from tests.context import pathtools


@ddt
@unittest.skip("slow test - internet connection required")
class FromEntrezFastaProviderTest_FastaCache(unittest.TestCase):

    def test_1_download_sequence_and_save_to_cache(self):
        fasta_provider = FromEntrezFastaProvider(email_address="paulina-ph@wp.pl", use_cache=True)
        cache_dir_path = pathtools.get_child_path(Path.cwd(), ".fastacache")
        if cache_dir_path.exists():
            shutil.rmtree(cache_dir_path)

        sequence_id = "AB050936.1"
        actual_sequence = fasta_provider.get_source(sequence_id, None, None)

        # cache directory creation
        cache_directory_created = cache_dir_path.exists()
        self.assertTrue(cache_directory_created)

        # file creation
        files_in_cache_dircetory = [*cache_dir_path.glob("*")]
        expected_filepath = pathtools.get_child_path(cache_dir_path, f"{sequence_id}.fasta")
        file_created_in_cache = expected_filepath in files_in_cache_dircetory
        self.assertTrue(file_created_in_cache)

        # file content
        control_fasta_path = pathtools.get_child_path(Path(os.getcwd()),
                                     "FastaProviders_Tests/FromEntrezFastaProvider_Tests/AB050936.1.fasta")
        with open(control_fasta_path) as fasta_file_hanlder:
            expected_content = fasta_file_hanlder.read()
        with open(expected_filepath) as fasta_file_handler:
            actual_content = fasta_file_handler.read()
        self.assertEqual(expected_content, actual_content)

    def test_2_read_seqeunce_from_cache_instead_downloading(self):
        fasta_provider = FromEntrezFastaProvider(email_address="paulina-ph@wp.pl", use_cache=True)
        cache_dir_path = pathtools.get_child_path(Path.cwd(), ".fastacache")
        if cache_dir_path.exists():
            shutil.rmtree(cache_dir_path)

        cache_dir_path.mkdir()
        sequence_id = "seq1"
        expected_sequence = "ctg"
        fake_fasta_path = pathtools.get_child_path(cache_dir_path, f"{sequence_id}.fasta")
        with open(fake_fasta_path, 'w') as fake_fasta_handler:
            fake_fasta_handler.write(f">{sequence_id} cached\n{expected_sequence}")

        actual_sequence = fasta_provider.get_source(sequence_id, None, None)

        self.assertEqual(expected_sequence, actual_sequence)


if __name__ == '__main__':
    unittest.main()
