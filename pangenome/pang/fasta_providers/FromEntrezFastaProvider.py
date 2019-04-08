import os
from pathlib import Path
from typing import NewType

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pangenome.pang.fasta_providers.FastaProvider import FastaProvider
from pangenome.pang.pangraph.custom_types import Sequence
from pangenome.pang.tools import loggingtools, pathtools

EntrezSequenceID = NewType("EntrezSequenceID", str)

detailed_logger = loggingtools.get_logger("details")


class FromEntrezFastaProvider(FastaProvider):
    def __init__(self, email_address: str, use_cache: bool):
        super().__init__()
        Entrez.email = email_address
        self.fasta_cache = FastaCache(Path(os.getcwd()))
        self.use_cache = use_cache

    def get_source(self, sequence_id: EntrezSequenceID, start: int = None, end: int = None) -> str:
        sequence_is_cached = self.fasta_cache.sequence_is_cached(sequence_id)
        if self.use_cache and sequence_is_cached:
            sequence = self.fasta_cache.read_from_cache(sequence_id)
        elif self.use_cache and not sequence_is_cached:
            sequence = self._download_from_ncbi(sequence_id, start, end)
            self.fasta_cache.save_to_cache(sequence_id, sequence)
        else:
            sequence = self._download_from_ncbi(sequence_id, start, end)
        return sequence

    def _download_from_ncbi(self, sequence_id: EntrezSequenceID, start: int, end: int) -> Sequence:
        detailed_logger.info(f"Downloading from entrez sequence {sequence_id}...")
        try:
            if start is not None and end is not None:
                handle = Entrez.efetch(db="nucleotide",
                                       id=sequence_id,
                                       rettype="fasta",
                                       retmode="text",
                                       seq_start=start,
                                       seq_stop=end)
            else:
                handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype="fasta", retmode="text")
            fasta_content = self.get_raw_sequence_from_fasta(handle)
            return fasta_content
        except Exception as ex:
            raise Exception(f"Cannot download from Entrez sequence of ID: {sequence_id}") from ex


class FastaCache:
    def __init__(self, parent_dir: Path):
        self.parent_dir = parent_dir
        self.cache_dir = pathtools.get_child_path(parent_dir, ".fastacache")

    def cache_dir_exists(self) -> bool:
        return pathtools.dir_exists(self.cache_dir)

    def create_cache_dir(self) -> None:
        if not self.cache_dir_exists():
            pathtools.create_dir(self.cache_dir)
            detailed_logger.info(".fastacache directory was created.")
        else:
            detailed_logger.warning("Cannot create .fastacache directory, as it already exists.")

    def save_to_cache(self, seq_id: EntrezSequenceID, sequence: Sequence)-> None:
        detailed_logger.info(f"Caching sequence {seq_id}...")
        if not self.cache_dir_exists():
            self.create_cache_dir()
        cache_filename = self.get_cached_filepath(seq_id)
        with open(cache_filename, 'w') as fasta_file_handle:
            SeqIO.write(SeqRecord(seq=Seq(sequence), id=seq_id, description="cached"), fasta_file_handle, "fasta")

    def read_from_cache(self, seq_id: EntrezSequenceID) -> Sequence:
        detailed_logger.info(f"Reading {seq_id} from cache...")
        cache_filepath = self.get_cached_filepath(seq_id)
        with open(cache_filepath) as fasta_handle:
            seq = SeqIO.read(fasta_handle, "fasta")
        return Sequence(seq.seq)

    def get_cached_filepath(self, seq_id: EntrezSequenceID) -> Path:
        return pathtools.get_child_path(self.cache_dir, f"{seq_id}.fasta")

    def sequence_is_cached(self, sequence_id: EntrezSequenceID) -> bool:
        if not self.cache_dir_exists():
            return False
        expected_fasta_file_name = self.get_cached_filepath(sequence_id)
        fasta_files_in_cache_dir = self.cache_dir.glob("*.fasta")
        if expected_fasta_file_name in [*fasta_files_in_cache_dir]:
            return True
        return False
