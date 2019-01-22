# -*- coding: utf-8 -*-
class SequenceInfo:
    def __init__(self, block_id, seq_id, start_pos, strand):
        self.block = block_id
        self.id = seq_id
        self.start = start_pos
        self.strand = strand