# -*- coding: utf-8 -*-
import six
from .graph.EdgeInfo import EdgeInfo
from collections import defaultdict, namedtuple

Edge = namedtuple('Edge', ['left', 'right', 'type'])

def edge_type(strand1, strand2):
    return (strand1, -strand2)
    
def connect_blocks(seq1, seq2, d):
    edgeType = edge_type(seq1.strand, seq2.strand)
    if seq1.block == seq2.block: 
        edgeType = edge_type(min(seq1.strand, seq2.strand), max(seq1.strand, seq2.strand))    
    elif seq1.block > seq2.block:
        edgeType = edgeType[::-1]
        seq1, seq2 = seq2, seq1   
    e = Edge(seq1.block, seq2.block, edgeType)
    d[e][0].append((EdgeInfo(seq1.id, seq1.start), EdgeInfo(seq2.id, seq2.start)))
    d[e][1] += 1

def weight(seq):
    d = defaultdict(lambda: [[], 0])
    if six.PY2:
        for i in xrange(len(seq) - 1):
            if seq[i].id == seq[i+1].id:
                connect_blocks(seq[i], seq[i+1], d)
    else:
        for i in range(len(seq) - 1):
            if seq[i].id == seq[i+1].id:
                connect_blocks(seq[i], seq[i+1], d)
    return d
