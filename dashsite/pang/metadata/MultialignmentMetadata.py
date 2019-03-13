import pandas

from pangraph.custom_types import SequenceID
from typing import Dict, List
from io import StringIO

class MultialignmentMetadata:
    metadata_df: Dict[SequenceID, Dict]

    def __init__(self,
                 metadata_file_content: str = None):
        if metadata_file_content is None:
            self.metadata_df = None
        else:
            self.metadata_df = self._read_metadata_csv(metadata_file_content)

    def _read_metadata_csv(self, csv_content):
        try:
            metadata_df = pandas.read_csv(StringIO(csv_content),sep=',',error_bad_lines=False)
        except Exception as e:
            raise Exception("Error when reading csv metadata.") from e

        try:
            metadata_df = metadata_df.set_index('seqid')
        except Exception as e:
            raise Exception("No \'seqid\' column in csv metadata.")

        return metadata_df

    def get_all_sequences_ids(self):
        return self.metadata_df.index.tolist()

    def feed_with_maf_data(self, names_in_maf: List[str]) -> None:
        self.metadata_df['mafname'] = self.metadata_df.index.map(lambda seqid: self._get_mafname(seqid, names_in_maf))
        pass

    def get_seq_metadata_as_dict(self, seq_id):
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
