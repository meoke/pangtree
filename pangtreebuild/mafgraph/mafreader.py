# -*- coding: utf-8 -*-
from Bio import AlignIO
from .graph.Block import Block
from .graph.SequenceInfo import SequenceInfo

def start_position(sequence):
    # Return start position relative to the plus strand
    if sequence.annotations["strand"] == -1:
        return sequence.annotations["srcSize"] - sequence.annotations["start"] - sequence.annotations["size"]
    else:
        return sequence.annotations["start"]
    
def read_maf(maf_file):
    blocks, seq = [], []
    for i, mafblock in enumerate(AlignIO.parse(maf_file, "maf")):
        blocks.append(Block(i, mafblock))
        for sequence in mafblock:
            seq.append(SequenceInfo(i, sequence.id, start_position(sequence), sequence.annotations["strand"]))
    seq.sort(key = lambda s: (s.id, s.start))
    return blocks, seq
