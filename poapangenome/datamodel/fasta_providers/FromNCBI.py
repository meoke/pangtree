import os
import re
from pathlib import Path
from typing import NewType, Dict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from datamodel.Node import Base
from datamodel.Sequence import SequenceID
from datamodel.fasta_providers.FastaProvider import FastaProvider, FastaProviderException
from Bio import Entrez, SeqIO

from tools import pathtools, logprocess
NCBISequenceID = NewType("NCBISequenceID", str)
detailed_logger = logprocess.get_logger("details")


class EmailAddress:
    """E-mail address requiered when Fasta Provider is \"NCBI\"
       as Entrez API obligates the user to pass e-mail address."""

    def __init__(self, email_address: str):
        match = re.match(r'^[_a-z0-9-]+(\.[_a-z0-9-]+)*@[a-z0-9-]+(\.[a-z0-9-]+)*(\.[a-z]{2,4})$', email_address)
        if match is not None:
            self.value = email_address
        else:
            raise FastaProviderException(f"Incorrect e-mail address ({email_address}).")


class FromNCBI(FastaProvider):
    def __init__(self, email_address: EmailAddress, use_cache: bool):
        Entrez.email = email_address.value
        self.fasta_disk_cache = FastaDiskCache(Path(os.getcwd()))
        self.use_cache = use_cache
        self.sequences: Dict[SequenceID, str] = {}

    def get_base(self, sequence_id: SequenceID, i: int) -> Base:
        if sequence_id not in self.sequences.keys():
            sequence_is_cached = self.fasta_disk_cache.sequence_is_cached(sequence_id)
            if self.use_cache and sequence_is_cached:
                self.sequences[sequence_id] = self.fasta_disk_cache.read_from_cache(sequence_id)
            elif self.use_cache and not sequence_is_cached:
                sequence = self._download_from_ncbi(sequence_id)
                self.sequences[sequence_id] = sequence
                self.fasta_disk_cache.save_to_cache(sequence_id, sequence)
            else:
                self.sequences[sequence_id] = self._download_from_ncbi(sequence_id)

        return Base(self.sequences[sequence_id][i])

    def _download_from_ncbi(self, sequence_id: SequenceID) -> str:
        detailed_logger.info(f"Downloading from entrez sequence {sequence_id}...")
        entrez_sequence_id = self._guess_ncbi_sequence_id(sequence_id)
        try:
            handle = Entrez.efetch(db="nucleotide", id=entrez_sequence_id, rettype="fasta", retmode="text")
            fasta_content = FastaProvider.get_raw_sequence_from_fasta(handle)
            return fasta_content
        except Exception as ex:
            raise Exception(f"Cannot download from Entrez sequence of ID: {sequence_id}") from ex

    def _guess_ncbi_sequence_id(self, seqid: SequenceID) -> str:
        detailed_logger.info(f"Guessing entrez sequence id...")
        version_indications = [*re.finditer('v[0-9]', seqid.value)]
        if len(version_indications) == 1:
            version_start = version_indications[0].span()[0]
            if version_start == len(seqid.value) - 2:
                guessed_entrez_name = seqid.value[0:version_start] + "." + seqid.value[version_start+1:]
            else:
                guessed_entrez_name = seqid.value
        else:
            guessed_entrez_name = seqid.value
        detailed_logger.info(f"{seqid} translated to {guessed_entrez_name}")
        return guessed_entrez_name

class FastaDiskCache:
    def __init__(self, parent_dir: Path):
        self.parent_dir = parent_dir
        self.cache_dir = pathtools.get_child_path(parent_dir, ".fastacache")

    def sequence_is_cached(self, sequence_id: SequenceID) -> bool:
        if not self.cache_dir_exists():
            return False
        expected_fasta_file_name = self.get_cached_filepath(sequence_id)
        fasta_files_in_cache_dir = self.cache_dir.glob("*.fasta")
        if expected_fasta_file_name in [*fasta_files_in_cache_dir]:
            return True
        return False

    def cache_dir_exists(self) -> bool:
        return pathtools.dir_exists(self.cache_dir)

    def get_cached_filepath(self, seq_id: SequenceID) -> Path:
        return pathtools.get_child_path(self.cache_dir, f"{seq_id}.fasta")

    def read_from_cache(self, seq_id: SequenceID) -> str:
        detailed_logger.info(f"Reading {seq_id} from cache...")
        cache_filepath = self.get_cached_filepath(seq_id)
        with open(cache_filepath) as fasta_handle:
            seq = SeqIO.read(fasta_handle, "fasta")
        return seq.seq

    def create_cache_dir(self) -> None:
        if not self.cache_dir_exists():
            pathtools.create_dir(self.cache_dir)
            detailed_logger.info(".fastacache directory was created.")
        else:
            detailed_logger.warning(".fastacache directory not created, as it already exists.")

    def save_to_cache(self, seq_id: SequenceID, sequence: str) -> None:
        detailed_logger.info(f"Caching sequence {seq_id}...")
        if not self.cache_dir_exists():
            self.create_cache_dir()
        cache_filename = self.get_cached_filepath(seq_id)
        with open(cache_filename, 'w') as fasta_file_handle:
            SeqIO.write(SeqRecord(seq=Seq(sequence), id=seq_id.value, description="cached"), fasta_file_handle, "fasta")


