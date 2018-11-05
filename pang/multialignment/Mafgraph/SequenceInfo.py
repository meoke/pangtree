class SequenceInfo:
    def __init__(self, sequence_id: int, start_pos: int):
        self.ID = sequence_id
        self.start_pos = start_pos

    def __str__(self):
        return f"Sequence {self.ID} starting at {self.start_pos}"
