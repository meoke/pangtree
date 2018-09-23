class SequenceMetadata:
    def __init__(self, genbankID: str, assemblyID: str, mafname: str, name: str, group: str):
        self.genbankID = genbankID
        self.assemblyID = assemblyID
        self.mafname = mafname
        self.name = name
        self.group = group
