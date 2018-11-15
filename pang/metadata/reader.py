import json
from .SequenceMetadata import SequenceMetadata
from .MultialignmentMetadata import MultialignmentMetadata


def read(gen_info_file, as_string=False):
    """Parse genomes metadata file to Metadata object."""

    if as_string == False:
        with open(gen_info_file) as p:
            metadata = json.load(p)
    else:
        metadata = json.loads(gen_info_file)
    genome_metadata = {seq['id']: SequenceMetadata(name=seq['name'],
                                                       genbankID=seq['genbankID'],
                                                       assemblyID=seq['assemblyID'],
                                                       mafname=seq['mafname'],
                                                       group=seq['group'])for seq in metadata['data']}
    return MultialignmentMetadata(title=metadata['title'],
                                      version=metadata['source'],
                                      genomes_metadata=genome_metadata)
