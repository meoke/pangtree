# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 07:59:21 2018

@author: Ania
"""

from Arc import Arc

class Block2:
    def __init__(self, block_id, alignment):
        self.id = block_id # nr bloku w pliku wejsciowym
        self.alignment = alignment 
        self.group = self # identyfikator zbioru, do ktorego nalezy blok
        self.ord = {self: 0} # slownik pozycji blokow nalezacych do tego samego zbioru 
        self.orient =  1  # orientacja danego bloku 
        self.out_edges = []
        self.max = 0 # najwieksza pozycja w order
        self.min = 0 # najmniejsza pozycja w order
        
    def orientation(self):
        return self.orient
    
    def find(self):
        return self.group
    
    def minimum(self):
        root = self.find()
        return root.min
    
    def maximum(self):
        root = self.find()
        return root.max
    
    def order(self):
        root = self.find()
        return root.ord[self]
    
    def size(self):
        root = self.find()
        return root.max - root.min + 1
    
    def reorder(self, i):
        root = self.find()
        root.ord[self] = i
    
    def unionto(self, other, reverse, flank):
        selfroot = self.find()
        otherroot = other.find()
        if reverse == -1:
            selfroot.max, selfroot.min = -selfroot.min, -selfroot.max
        if flank == -1:
            n = otherroot.min - selfroot.max - 1
            otherroot.min -= selfroot.size()
        else:
            n = otherroot.max - selfroot.min + 1
            otherroot.max += selfroot.size()   
        for block in selfroot.ord:
            selfroot.ord[block] *= reverse
            block.orient *= reverse
            block.group = otherroot
            otherroot.ord[block] = n + selfroot.ord[block]
        del selfroot.ord
        del selfroot.max
        del selfroot.min     
       
    def add_out_edges(self, to, list_of_seq):
        self.out_edges.append(Arc(to, list_of_seq))
