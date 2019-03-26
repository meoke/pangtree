import re

import pandas as pd

from fasta_providers.FromEntrezFastaProvider import EntrezSequenceID
from pangraph.custom_types import SequenceID
from typing import Dict, List, Any
from io import StringIO

from tools import loggingtools

global_logger = loggingtools.get_global_logger()
detailed_logger = loggingtools.get_logger("details")


class MultialignmentMetadata:
    metadata_df: pd.DataFrame

    def __init__(self,
                 metadata_file_content: str = None):
        if metadata_file_content is None:
            self.metadata_df = None
        else:
            self.metadata_df = self._read_metadata_csv(metadata_file_content)

    def _read_metadata_csv(self, csv_content):
        global_logger.info("Reading metadata...")

        csv_erros = self.check_csv_correctness(csv_content)
        if csv_erros is not None:
            raise Exception(csv_erros)

        try:
            metadata_df = pd.read_csv(StringIO(csv_content), sep=',', error_bad_lines=True)
        except Exception as e:
            raise Exception("Error when reading csv metadata.") from e

        try:
            metadata_df = metadata_df.set_index('seqid')
        except Exception as e:
            raise Exception("No \'seqid\' column in csv metadata.")

        seqid_column_values = metadata_df.index
        if sorted(set(seqid_column_values)) != sorted(seqid_column_values):
            raise Exception("Not unique values seqid column in metadata file. Make them unique.")

        return metadata_df

    def get_all_sequences_ids(self):
        return [SequenceID(seq_id) for seq_id in self.metadata_df.index.tolist()]

    def feed_with_maf_data(self, names_in_maf: List[str]) -> None:
        seq_id_to_full_mafname = [{'seqid': MultialignmentMetadata.get_seqid_from_mafname(mafname), 'mafname': mafname}
                                  for mafname in names_in_maf]
        if self.metadata_df is None:
            self.metadata_df = pd.DataFrame.from_dict(seq_id_to_full_mafname)
            self.metadata_df = self.metadata_df.set_index('seqid')
        else:
            if 'mafname' not in list(self.metadata_df):
                self.metadata_df['mafname'] = None
            for seq_id_mafname in seq_id_to_full_mafname:
                seq_id = seq_id_mafname['seqid']
                mafname = seq_id_mafname['mafname']
                if seq_id in self.metadata_df.index:
                    self.metadata_df.loc[seq_id, 'mafname'] = mafname
                else:
                    self.metadata_df.loc[seq_id, 'mafname'] = mafname
        # todo check if no column contains "CONSENSUS" - do we need this?

    def get_seq_metadata_as_dict(self, seq_id: SequenceID) -> Dict[str, Any]:
        try:
            return self.metadata_df.to_dict(orient='index')[seq_id]
        except Exception as e:
            raise Exception(f"No record for {seq_id} in metadata.") from e

    def _get_mafname(self, seqid: SequenceID, names_in_maf: List[str]) -> str:
        for n in names_in_maf:
            splitted = n.split('.')
            if (len(splitted) > 1 and splitted[1] == seqid) or splitted[0] == seqid:
                return n
        return None

    @staticmethod
    def get_seqid_from_mafname(mafname):
        splitted = mafname.split('.')
        if len(splitted) > 1:
            return ".".join(splitted[1:])
        elif len(splitted) == 1:
            return splitted[0]

    def get_entrez_name(self, seqid: SequenceID) -> EntrezSequenceID:
        try:
            return EntrezSequenceID(self.metadata_df[seqid]['entrez'])
        except KeyError:
            return self._guess_entrez_name(seqid)

    def _guess_entrez_name(self, seqid: SequenceID) -> EntrezSequenceID:
        detailed_logger.info(f"Guessing entrez sequence id...")
        version_indications = [*re.finditer('v[0-9]{1}', seqid)]
        guessed_entrez_name = seqid
        if len(version_indications) == 1:
            version_start = version_indications[0].span()[0]
            if version_start == len(seqid) - 2:
                guessed_entrez_name = seqid[0:version_start] + "." + seqid[version_start+1:]
        detailed_logger.info(f"{seqid} translated to {guessed_entrez_name}")
        return EntrezSequenceID(guessed_entrez_name)

    def check_csv_correctness(self, csv_content):
        if not csv_content:
            return "Empty csv file."
        csv_content = StringIO(csv_content)
        header = csv_content.readline().split(',')
        for i, line in enumerate(csv_content):
            if len(line.split(',')) != len(header):
                return f"CSV metadata error. Different fields number in line {i} than in header line."
        return None
