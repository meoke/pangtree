class SequenceMetadata:
    def __init__(self,
                 name: str,
                 genbankID: str = "",
                 assemblyID: str = "",
                 mafname: str = "",
                 group: str = ""):
        self.name = name
        self.genbankID = genbankID
        self.assemblyID = assemblyID
        self.mafname = mafname
        self.group = group
