import json
from .SequenceMetadata import SequenceMetadata
from .MultialignmentMetadata import MultialignmentMetadata


def read(metadata_json: str) -> MultialignmentMetadata:
    """Parse genomes metadata json to MultialignmentMetadata object."""

    metadata_json = json.loads(metadata_json)
    genome_metadata = {seq['mafname']: SequenceMetadata(name=seq['name'],
                                                        genbankID=seq['genbankID'],
                                                        assemblyID=seq['assemblyID'],
                                                        mafname=seq['mafname'],
                                                        group=seq['group'])
                       for seq in metadata_json['data']}

    return MultialignmentMetadata(title=metadata_json['title'],
                                  version=metadata_json['source'],
                                  genomes_metadata=genome_metadata)

