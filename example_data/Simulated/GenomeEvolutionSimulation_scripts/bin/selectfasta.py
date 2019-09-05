#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 17:29:25 2019

@author: norbert
"""

from Bio import SeqIO
from sys import argv


seqtype = argv[1]
seqnames = argv[2]
dirname = argv[3]

#specnames = 'ABCDEFGHIJ'
#chrname = 'chr20'


records = []

seqlist = seqnames.split(',')
#print seqlist
for seqname in seqlist:
#    print seqname
    s,c=seqname.split('.')
    for record in SeqIO.parse(dirname+s+'/seq.fa', "fasta"):
        if record.id==c:
            record.id=seqname
            record.description = ''
            record.name = ''
            records.append(record)
SeqIO.write(records,dirname+seqtype+".fasta","fasta")

