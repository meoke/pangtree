import os
import shutil
import unittest
from pathlib import Path

from ddt import ddt

from tests.context import FromNCBI
from tests.context import pathtools
from tests.context import pSeq, pNode


@ddt
class FromNCBI_FastaDiskCache_Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_provider = FromNCBI(use_cache=False)

    @unittest.skip('Internet connection required -> long execution')
    def test_1_download_sequence_and_save_to_cache(self):
        fasta_provider = FromNCBI(use_cache=True)
        cache_dir_path = pathtools.get_child_path(Path.cwd(), ".fastacache")
        if cache_dir_path.exists():
            shutil.rmtree(cache_dir_path)

        sequence_id = pSeq.SequenceID("AB050936.1", skip_part_before_dot=False)
        actual_sequence = fasta_provider.get_base(sequence_id, 15)

        # cache directory creation
        cache_directory_created = cache_dir_path.exists()
        self.assertTrue(cache_directory_created)

        # file creation
        files_in_cache_dircetory = [*cache_dir_path.glob("*")]
        expected_filepath = pathtools.get_child_path(cache_dir_path, f"{sequence_id}.fasta")
        file_created_in_cache = expected_filepath in files_in_cache_dircetory
        self.assertTrue(file_created_in_cache)

        # file content
        control_fasta_path = Path('tests/data/fasta_providers/fasta_ncbi/AB050936.1.fasta')

        with open(control_fasta_path) as fasta_file_hanlder:
            expected_content = fasta_file_hanlder.read()
        with open(expected_filepath) as fasta_file_handler:
            actual_content = fasta_file_handler.read()
        self.assertEqual(expected_content, actual_content)

    def test_2_read_seqeunce_from_cache_instead_downloading(self):
        fasta_provider = FromNCBI(use_cache=True)
        cache_dir_path = pathtools.get_child_path(Path.cwd(), ".fastacache")
        if cache_dir_path.exists():
            shutil.rmtree(cache_dir_path)

        cache_dir_path.mkdir()
        sequence_id = pSeq.SequenceID("seq1")
        fake_sequence = "foo"
        expected_base = pNode.Base("o")
        fake_fasta_path = pathtools.get_child_path(cache_dir_path, f"{sequence_id}.fasta")
        with open(fake_fasta_path, 'w') as fake_fasta_handler:
            fake_fasta_handler.write(f">{sequence_id} cached\n{fake_sequence}")

        actual_base = fasta_provider.get_base(sequence_id, 2)

        self.assertEqual(expected_base, actual_base)


if __name__ == '__main__':
    unittest.main()
