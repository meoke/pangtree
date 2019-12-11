import shutil
import unittest
from pathlib import Path

from ddt import ddt

from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import missings, msa
from pangtreebuild.tools import pathtools


@ddt
class FromNCBI_FastaDiskCache_Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.fasta_provider = missings.FromNCBI(use_cache=False)

    # @unittest.skip('Internet connection required -> long execution')
    def test_1_download_sequence_and_save_to_cache(self):
        cache_dir_path = pathtools.get_child_path(Path.cwd(), ".fastacache")
        if cache_dir_path.exists():
            shutil.rmtree(cache_dir_path)

        ncbi_fasta_provider = missings.FromNCBI(use_cache=True)
        sequence_id = msa.SequenceID("AB050936v1")

        _ = ncbi_fasta_provider.get_base(sequence_id, 0)

        # cache directory creation
        cache_directory_created = cache_dir_path.exists()
        self.assertTrue(cache_directory_created)

        # file creation
        files_in_cache_dircetory = [*cache_dir_path.glob("*")]
        expected_filepath = pathtools.get_child_path(cache_dir_path, f"{sequence_id}.fasta")
        file_created_in_cache = expected_filepath in files_in_cache_dircetory
        self.assertTrue(file_created_in_cache)

        # file content
        control_fasta_path = Path(__file__).parent.joinpath('fasta_ncbi/AB050936.1.fasta').resolve()

        with open(control_fasta_path) as fasta_file_hanlder:
            expected_content = fasta_file_hanlder.read()
        with open(expected_filepath) as fasta_file_handler:
            actual_content = fasta_file_handler.read()
        self.assertEqual(expected_content, actual_content)

    def test_2_read_seqeunce_from_cache_instead_downloading(self):
        fasta_provider = missings.FromNCBI(use_cache=True)
        cache_dir_path = pathtools.get_child_path(Path.cwd(), ".fastacache")
        if cache_dir_path.exists():
            shutil.rmtree(cache_dir_path)

        cache_dir_path.mkdir()
        sequence_id = msa.SequenceID("seq1")
        fake_sequence = "foo"
        expected_base = graph.Base("o")
        fake_fasta_path = pathtools.get_child_path(cache_dir_path, f"{sequence_id}.fasta")
        with open(fake_fasta_path, 'w') as fake_fasta_handler:
            fake_fasta_handler.write(f">{sequence_id} cached\n{fake_sequence}")

        actual_base = fasta_provider.get_base(sequence_id, 2)

        self.assertEqual(expected_base, actual_base)


if __name__ == '__main__':
    unittest.main()
