import os
from pathlib import Path
from typing import NewType

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fasta_providers.FastaProvider import FastaProvider
from tools import loggingtools, pathtools

EntrezSequenceID = NewType("EntrezSequenceID", str)

detailed_logger = loggingtools.get_logger("details")


class FromEntrezFastaProvider(FastaProvider):
    def __init__(self, email_address: str, use_cache: bool):
        super().__init__()
        Entrez.email = email_address
        self.fasta_cache = FastaCache(Path(os.getcwd()))
        self.use_cache = use_cache

    def get_source(self, sequenceID: EntrezSequenceID, start: int = None, end: int = None) -> str:
        sequence_is_cached = self.fasta_cache.sequence_cached(sequenceID)
        if self.use_cache and sequence_is_cached:
            sequence = self.fasta_cache.read_from_cache(sequenceID)
        elif self.use_cache and not sequence_is_cached:
            sequence = self.download_from_ncbi(sequenceID, start, end)
            self.fasta_cache.save_to_cache(sequenceID, sequence)
        else:
            sequence = self.download_from_ncbi(sequenceID, start, end)
        return sequence

    def download_from_ncbi(self, sequenceID, start, end):
        detailed_logger.info(f"Downloading from entrez sequence {sequenceID}...")
        try:
            if start is not None and end is not None:
                handle = Entrez.efetch(db="nucleotide",
                                       id=sequenceID,
                                       rettype="fasta",
                                       retmode="text",
                                       seq_start=start,
                                       seq_stop=end)
            else:
                handle = Entrez.efetch(db="nucleotide", id=sequenceID, rettype="fasta", retmode="text")
            fasta_content = self.get_raw_sequence_from_fasta(handle)
            return fasta_content
        except Exception as ex:
            raise Exception(f"Cannot download from Entrez sequence of ID: {sequenceID}") from ex


class FastaCache:
    def __init__(self, parent_dir: Path):
        self.parent_dir = parent_dir
        self.cache_dir = pathtools.get_child_path(parent_dir, ".fastacache")

    def cache_dir_exists(self):
        return pathtools.dir_exists(self.cache_dir)

    def create_cache_dir(self) -> None:
        if not self.cache_dir_exists():
            pathtools.create_dir(self.cache_dir)
            detailed_logger.info(".fastacache directory was created.")
        else:
            detailed_logger.warning("Cannot create .fastacache directory, as it already exists.")

    def save_to_cache(self, seq_id, sequence)-> None:
        detailed_logger.info(f"Caching sequence {seq_id}...")
        cache_filename = self.get_cached_filepath(seq_id)
        with open(cache_filename, 'w') as fasta_file_handle:
            SeqIO.write(SeqRecord(seq=Seq(sequence), id=seq_id, description="cached"), fasta_file_handle, "fasta")

    def read_from_cache(self, seq_id):
        detailed_logger.info(f"Reading {seq_id} from cache...")
        cache_filename = self.get_cached_filepath(seq_id)
        seq = SeqIO.read(cache_filename, "fasta")
        return seq.seq

    def get_cached_filepath(self, seq_id):
        return pathtools.get_child_path(self.cache_dir, self.get_cached_filename(seq_id))

    def get_cached_filename(self, seq_id):
        return f"{seq_id}.fasta"

    def sequence_cached(self, sequenceID):
        if not self.cache_dir_exists():
            return False
        expected_fasta_file_name = self.get_cached_filename(sequenceID)
        fasta_files_in_cache_dir = self.cache_dir.glob("*.fasta")
        if expected_fasta_file_name in fasta_files_in_cache_dir:
            return True
        return False