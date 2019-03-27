class SequenceMetadata:
    def __init__(self,
                 name: str,
                 genbank_id: str = "",
                 assembly_id: str = "",
                 mafname: str = "",
                 group: str = ""):
        self.name = name
        self.genbank_id = genbank_id
        self.assembly_id = assembly_id
        self.mafname = mafname
        self.group = group
