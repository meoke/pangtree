"""Missing in MAF nucleotide or protein providers."""

import abc
from io import StringIO
import os
from pathlib import Path
import re
from typing import Dict, Optional
import zipfile

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from pangtreebuild.pangenome import graph
from pangtreebuild.pangenome.parameters import msa
from pangtreebuild.tools import pathtools, logprocess

detailed_logger = logprocess.get_logger("details")


class FastaProviderException(Exception):
    """Any exception connected with providing FASTA sequences."""

    pass


class FastaProvider(abc.ABC):
    """Base class for FASTA providers.

    MAF file usually contains not full sequences but only the aliggned parts.
    To build an exact poagraph full sequences must be retrieved
    from NCBI or FASTA file."""

    @abc.abstractmethod
    def get_base(self, seq_id: msa.SequenceID, i: int) -> graph.Base:
        """Returns base at position i in sequence identified by seq_id.

        Args:
            seq_id: Sequence ID of the sequence being currently complemented.
            i: position of the base in the sequence.

        Returns:
            The nucleotide or protein base.
        """

        pass

    @staticmethod
    def get_sequence_from_fasta(fasta_content: StringIO) -> str:
        """Returns sequence from FASTA file.

        FASTA file contains a header line and the sequence
        splitted into several lines. This method returns just the sequence.

        Args:
            fasta_content: Content of FASTA file.

        Returns:
            Return sequence saved in fasta structured file.
        """

        _ = fasta_content.readline()
        return fasta_content.read().replace('\n', '')


class MissingBase(object):
    """A constant symbol (default "?") replacing sequence missing symbols.

    Args:
        symbol: one-character long str to be saved as the missing symbol

    Attributes:
        value (str): the symbol as str
    """

    def __init__(self, symbol: Optional[str] = '?'):
        if len(symbol) != 1:
            raise ValueError('Missing symbol must be a single character.')
        self.value = symbol


class ConstBaseProvider(FastaProvider):
    """Missing base provider which returns always the same single character.

    Args:
        missing_base: the character to be always returned as missing base

    Attributes:
        missing_base (graph.Base): a constant sequence base
    """

    def __init__(self, missing_base: MissingBase):
        self.missing_base: graph.Base = graph.Base(missing_base.value)

    def get_base(self, seq_id: msa.SequenceID, i: int) -> graph.Base:
        """Returns const base symbol.

        Args:
            seq_id: ignored as the returned symbol is always the same.
            i: ignored as the returned symbol is always the same.

        Returns:
            The const symbol serving as missing symbol.
        """

        return self.missing_base


class FromFile(FastaProvider):
    """Missing base provider with FASTA file or zipped FASTA files source.

    Args:
        fasta: Path to the FASTA file or zipped FASTA files.
    """

    def __init__(self, fasta: Path):
        self._sequences: Dict[msa.SequenceID, str] = self._read_fastas(fasta)

    def get_base(self, seq_id: msa.SequenceID, i: int) -> graph.Base:
        """Returns base at position i in sequence identified by seq_id.

        Args:
            seq_id: ID of sequence present in given fasta file.
            i: position of the base to check in fasta file.

        Returns:
            The base present in seuqence seq_id at position i.
        """

        if seq_id not in self._sequences.keys():
            raise FastaProviderException(f"Wrong sequence id: {seq_id}. ")
        if i > len(self._sequences[seq_id]):
            raise FastaProviderException(f"""Index {i} is too large
                                          for sequence {seq_id}.""")
        return graph.Base(self._sequences[seq_id][i])

    def _read_fastas(self, fastas_file: Path) -> Dict[msa.SequenceID, str]:
        """Reads FASTA content to internal dictionary.

        Args:
            fastas_file: Path to the fasta file or zipped fasta files.

        Returns:
            Dictionary of sequence_id to sequence content.

        Raises:
            FastaProviderException: If there is no extension or the extension
            is unknown (other than zip/fasta/fna/faa.
        """

        fastas_file_suffix = fastas_file.suffix
        if fastas_file_suffix == '':
            raise FastaProviderException("""No file extension in fasta source.
                                        Append extension to file name.
                                        Supported: zip/fasta/fna/faa.""")
        else:
            fastas_file_extensions = fastas_file_suffix[1:]
            if fastas_file_extensions == "zip":
                return self._read_zip(fastas_file)
            elif fastas_file_extensions in ["fna", "fasta", "faa"]:
                return self._read_fasta(fastas_file)
            else:
                raise FastaProviderException("""Unknown fasta file format.
                                             Supported: zip/fasta/fna/faa.""")

    def _read_zip(self, zip_path: Path) -> Dict[msa.SequenceID, str]:
        """For given zip of fastas reads the sequences into a dictionary.

        Args:
            zip_path: Path tp the zip.

        Returns:
            Dictionary of sequence_id to sequence content.

        Raises:
            FastaProviderException: If no seuquences were found or
                                    the zipped fasta files were incorrect.

        """

        if not zipfile.is_zipfile(zip_path):
            raise FastaProviderException("Incorrect zip fasta source.")
        fastas_dict = {}
        with zipfile.ZipFile(zip_path) as myzip:
            zipinfolist = myzip.infolist()
            for zipinfo in zipinfolist:
                with myzip.open(zipinfo) as myfile:
                    fasta = myfile.read().decode("ASCII")
                    for record in SeqIO.parse(StringIO(fasta), "fasta"):
                        self._add_record_to_dict(record, fastas_dict)
        if len(fastas_dict) == 0:
            raise FastaProviderException("No sequences in zipped fastas or incorrect zipped files.")
        return fastas_dict

    def _read_fasta(self, fasta_path: Path) -> Dict[msa.SequenceID, str]:
        """Read sequences content from a single fasta file.

        Args:
            fasta_path: Path to the fasta file.

        Returns:
            Dictionary of sequence_id to sequence content.

        Raises:
            FastaProviderException: If no sequences were found or
                                    the zipped fasta files were incorrect.
        """
        fastas_dict = {}
        with open(fasta_path) as fasta_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                self._add_record_to_dict(record, fastas_dict)
        if len(fastas_dict) == 0:
            raise FastaProviderException("No sequences in zipped fastas or incorrect zipped files.")
        return fastas_dict

    def _add_record_to_dict(self,
                            record: SeqRecord,
                            fastas_dict: Dict[msa.SequenceID, str]) -> None:
        """Adds record the sequence to dict if the sequence is not empty.

        Args:
            record: The sequence to add.

        Raises:
            FastaProviderExecption: If record contains empty sequence or
                                    its ID is already present in the dict.
        """
        if len(record.seq) == 0:
            raise FastaProviderException("Empty sequence in FASTA. Provide the sequence or remove its header.")
        if msa.SequenceID(str(record.id)) in fastas_dict.keys():
            raise FastaProviderException("Incorrect fasta provided: sequences IDs are not unique.")
        fastas_dict[msa.SequenceID(str(record.id))] = record.seq


class FromNCBI(FastaProvider):
    """Missing base provider with NCBI as source.

    Args:
        use_cache:  Switch whether to save locally sequences downloaded
                    from NCBI and reuse them in further calls.
    """

    def __init__(self, use_cache: bool):
        Entrez.email = "paulinahyzy@gmail.com"
        self._fasta_disk_cache = _FastaDiskCache(Path(os.getcwd()))
        self._use_cache = use_cache
        self._sequences: Dict[msa.SequenceID, str] = {}

    def get_base(self, seq_id: msa.SequenceID, i: int) -> graph.Base:
        """Returns base at position i in sequence identified in NCBI
        by seq_id or sth similar.

        Args:
            seq_id: sequence_id, used as is or some guessing is performed
                    if no results available.
            i: position of the base to check in fasta file.

        Returns:
            The base present in seuqence sequence_id at position i.
        """
        if seq_id not in self._sequences.keys():
            sequence_is_cached = self._fasta_disk_cache.seq_is_cached(seq_id)
            if self._use_cache and sequence_is_cached:
                self._sequences[seq_id] = self._fasta_disk_cache.read(seq_id)
            elif self._use_cache and not sequence_is_cached:
                sequence = self._download_from_ncbi(seq_id)
                self._sequences[seq_id] = sequence
                self._fasta_disk_cache._save_to_cache(seq_id, sequence)
            else:
                self._sequences[seq_id] = self._download_from_ncbi(seq_id)

        return graph.Base(self._sequences[seq_id][i])

    def _download_from_ncbi(self, seq_id: msa.SequenceID) -> str:
        """Connects to NCBI and downloads full sequence identified by seq_id.

        Args:
            seq_id: Sequence ID in format [seq_id] or [seq_id]v[version].

        Returns:
            Sequence as str.

        Raises:
            ConnectionError: If any error while connecting NCBI appeared.
        """

        detailed_logger.info(f"Downloading from entrez sequence {seq_id}...")
        entrez_sequence_id = self._guess_ncbi_sequence_id(seq_id)
        try:
            handle = Entrez.efetch(db="nucleotide",
                                   id=entrez_sequence_id,
                                   rettype="fasta",
                                   retmode="text")
            fasta_content = FastaProvider.get_sequence_from_fasta(handle)
            return fasta_content
        except Exception as ex:
            raise ConnectionError(f"""Cannot download from Entrez
                                    sequence of ID: {seq_id}""") from ex

    def _guess_ncbi_sequence_id(self, seqid: msa.SequenceID) -> str:
        """Tries to guess sequence id by searching for version indicator.

        Args:
            seqid: Sequence ID in format [seq_id] or [seq_id]v[version].

        Returns:
            Guessed sequence ID.
        """

        detailed_logger.info(f"Guessing entrez sequence id...")
        version_indications = [*re.finditer('v[0-9]', seqid.value)]
        if len(version_indications) == 1:
            v_start = version_indications[0].span()[0]
            if v_start == len(seqid.value) - 2:
                guessed_entrez_name = ".".join([seqid.value[0:v_start],
                                                seqid.value[v_start+1:]])
            else:
                guessed_entrez_name = seqid.value
        else:
            guessed_entrez_name = seqid.value
        detailed_logger.info(f"{seqid} translated to {guessed_entrez_name}")
        return guessed_entrez_name


class _FastaDiskCache(object):
    """Manages caching sequences downloaded from NCBI.

    Args:
        parent_dir: Path to the directory where new directory
                    with cached sequences will be created.
    """

    def __init__(self, parent_dir: Path):
        self._parent_dir = parent_dir
        self._cache_dir = pathtools.get_child_path(parent_dir, ".fastacache")

    def seq_is_cached(self, seq_id: msa.SequenceID) -> bool:
        """Returns information whether sequence with seq_id is cached.

        Args:
            seq_id: Sequence ID.

        Returns:
            True if sequence is cached, False if not.
        """

        if not self._cache_dir_exists():
            return False
        expected_fasta_file_name = self._get_cached_filepath(seq_id)
        fasta_files_in_cache_dir = self._cache_dir.glob("*.fasta")
        if expected_fasta_file_name in [*fasta_files_in_cache_dir]:
            return True
        return False

    def _cache_dir_exists(self) -> bool:
        return pathtools.dir_exists(self._cache_dir)

    def _get_cached_filepath(self, seq_id: msa.SequenceID) -> Path:
        return pathtools.get_child_path(self._cache_dir, f"{seq_id}.fasta")

    def read(self, seq_id: msa.SequenceID) -> str:
        """Gets full sequence identified by seq_id from cached local file.

        Args:
            seq_id: Sequence ID.

        Returns:
            The sequence as str.
        """

        detailed_logger.info(f"Reading {seq_id} from cache...")
        cache_filepath = self._get_cached_filepath(seq_id)
        with open(cache_filepath) as fasta_handle:
            seq = SeqIO.read(fasta_handle, "fasta")
        return seq.seq

    def _create_cache_dir(self) -> None:
        if not self._cache_dir_exists():
            pathtools.create_dir(self._cache_dir)
            detailed_logger.info(".fastacache directory was created.")
        else:
            detailed_logger.warning(""".fastacache directory not created,
                                    as it already exists.""")

    def _save_to_cache(self, seq_id: msa.SequenceID, sequence: str) -> None:
        detailed_logger.info(f"Caching sequence {seq_id}...")
        if not self._cache_dir_exists():
            self._create_cache_dir()
        cache_filename = self._get_cached_filepath(seq_id)
        with open(cache_filename, 'w') as fasta_file_handle:
            SeqIO.write(SeqRecord(seq=Seq(sequence),
                                  id=seq_id.value,
                                  description="cached"),
                        fasta_file_handle,
                        "fasta")
