# -*- coding: utf-8 -*-
import sys
sys.path.append('..')
from ..mafreader import start_position

class EdgeInfo:
    def __init__(self, seq_id, start_pos):
        self.seq_id = seq_id
        self.start = start_pos
        
    def set_start_position(self, alignment, orientation):
        for seq in alignment:
            if seq.id == self.seq_id and start_position(seq) == self.start:
                if orientation*seq.annotations["strand"] == 1:
                    break
                elif orientation == -1:
                    self.start = seq.annotations["srcSize"] - seq.annotations["start"] - seq.annotations["size"]
                    break
                else:
                    self.start = seq.annotations["start"]
                    break
