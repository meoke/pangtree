# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 07:00:31 2018

@author: Ania
"""
import networkx as nx
from Bio import AlignIO
from block2 import Block2
from collections import defaultdict, namedtuple

class SequenceInfo:
    def __init__(self, block_id, seq_id, start_pos, strand):
        self.block = block_id
        self.id = seq_id
        self.start = start_pos
        self.strand = strand

class EdgeInfo:
    def __init__(self, seq_id, start_pos):
        self.seq_id = seq_id
        self.start = start_pos

        
Edge = namedtuple('Edge', ['left', 'right', 'type'])
 
def start_position(sequence):
    # Return start position relative to the plus strand
    if sequence.annotations["strand"] == -1:
        return sequence.annotations["srcSize"] - sequence.annotations["start"] - sequence.annotations["size"]
    else:
        return sequence.annotations["start"]
    
def read_maf(maf_file):
    blocks, seq = [], []
    for i, mafblock in enumerate(AlignIO.parse(maf_file, "maf")):
        blocks.append(Block2(i, mafblock))
        for sequence in mafblock:
            seq.append(SequenceInfo(i, sequence.id, start_position(sequence), sequence.annotations["strand"]))
    seq.sort(key = lambda s: (s.id, s.start))
    return blocks, seq

def edge_type(strand1, strand2):
    return (strand1, -strand2)
    
def change_edge_type(edgeType):
    return edgeType[::-1]
   
def connect_blocks(seq1, seq2, d):
    edgeType = edge_type(seq1.strand, seq2.strand)
    if seq1.block == seq2.block: 
        edgeType = edge_type(min(seq1.strand, seq2.strand), max(seq1.strand, seq2.strand))    
    elif seq1.block > seq2.block:
        edgeType = change_edge_type(edgeType)
        seq1, seq2 = seq2, seq1   
    e = Edge(seq1.block, seq2.block, edgeType)
    d[e][0].append((EdgeInfo(seq1.id, seq1.start), EdgeInfo(seq2.id, seq2.start)))
    d[e][1] += 1   
    
def weight(seq):
    d = defaultdict(lambda: [[],0])
    # for i in xrange(len(seq)-1):
    for i in range(len(seq)-1):
        if seq[i].id == seq[i+1].id:
            connect_blocks(seq[i], seq[i+1], d)
    return d

""" TOPOLOGICAL SORTING """
def out_edges(v, blocks, ub, G):
    # Return a list of nodes connected to node v by edge leaving the vertex v   
    return [u for u in G[v] if ub >= blocks[u].order() > blocks[v].order()]

def in_edges(v, blocks, lb, G):
    # Return a list of nodes connected to node v by edge coming in to the vertex v     
    return [u for u in G[v] if lb < blocks[u].order() < blocks[v].order()]  

def dfs_f(root, ub, G, blocks):
    visited = set()
    stack = [root,]
    while stack:
        node = stack.pop()
        if blocks[node].order()==ub: return []
        if node not in visited:
            visited.add(node)
            stack.extend([x for x in out_edges(node, blocks, ub, G) if x not in visited])
    return list(visited) 

def dfs_b(root, lb, G, blocks):
    visited = set()
    stack = [root,]
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            stack.extend([x for x in in_edges(node, blocks, lb, G) if x not in visited])
    return list(visited)

def reorder(R_f, R_b, blocks):
    R_f.sort(key = lambda x: blocks[x].order()) 
    R_b.sort(key = lambda x: blocks[x].order())
    L = R_b + R_f
    O = sorted(blocks[x].order() for x in L)
    # for i in xrange(len(L)):
    for i in range(len(L)):
        blocks[L[i]].reorder(O[i])

def add_edge_within_component(x, y, G, blocks):
    lb = blocks[y].order()
    ub = blocks[x].order()
    if lb is ub: return
    elif lb < ub:
        R_f = dfs_f(y, ub, G, blocks)
        if R_f:
            R_b = dfs_b(x, lb, G, blocks)
            reorder(R_f, R_b, blocks)
            G.add_edge(x,y)
    else: 
        G.add_edge(x,y)
"""""" 

       
def add_edge_between_components(e, blocks):
    # Add edge between blocks of different components
    reverse, flank = 1, 1
    if blocks[e.left].size() < blocks[e.right].size():
        e = Edge(e.right, e.left, e.type[::-1])
    if blocks[e.left].orientation()*blocks[e.right].orientation() is e.type[0]*e.type[1]:
        reverse = -1
    if blocks[e.left].orientation()*e.type[0] < 0: 
        flank = -1
    blocks[e.right].unionto(blocks[e.left], reverse, flank)

def orient_block(block):
    if block.orientation() == -1:
        for u in block.alignment:
            u.seq = u.seq.reverse_complement()
            u.annotations["strand"] *= -1
            u.annotations["start"] = u.annotations["srcSize"] - u.annotations["size"] - u.annotations["start"]

def connect_components(blocks):
    d = {blocks[0].find(): 0}
    n = blocks[0].maximum()
    for block in blocks:
        if block.find() not in d:
            i = n - block.minimum() + 1
            d[block.find()] = i 
            block.reorder(i + block.order())
            n += block.size()
        else:
            block.reorder(d[block.find()] + block.order())
            
def sort_mafblocks(maf_file):   
    blocks, seq = read_maf(maf_file) # blocks - list of Block instances
    d = weight(seq)  # d - dictionary of edges
    edges = sorted(d.keys(), key = lambda x: d[x][1]) # list of edges sorted by the weight
    G = nx.Graph()
    for e in edges:
        if blocks[e.left].find() is blocks[e.right].find():
            if blocks[e.left].orientation()*blocks[e.right].orientation() is e.type[0]*e.type[1]:
                continue
            elif blocks[e.left].orientation()*e.type[0] > 0:
                add_edge_within_component(e.left, e.right, G, blocks)
            else:
                add_edge_within_component(e.right, e.left, G, blocks)           
        else:
            G.add_edge(e.left, e.right)
            add_edge_between_components(e, blocks)
    if len(blocks) != blocks[0].size(): 
        connect_components(blocks)
    else: 
        for block in blocks:
            orient_block(block)
    for edge in d:
        if (blocks[edge.left].order() < blocks[edge.right].order() 
         and (blocks[edge.left].orientation()*edge.type[0], blocks[edge.right].orientation()*edge.type[1]) == (1,-1)):
            blocks[edge.left].add_out_edges(edge.right, d[edge][0])
        elif (blocks[edge.left].order() > blocks[edge.right].order() 
         and (blocks[edge.right].orientation()*edge.type[1], blocks[edge.left].orientation()*edge.type[0]) == (1,-1)):
            sequences = [x[::-1] for x in d[edge][0]]
            blocks[edge.right].add_out_edges(edge.left, sequences)
    blocks = sorted(blocks, key = lambda b: b.order())
    return blocks


""" W powyższej pętli (for edge in d:) każdemu blokowi dodaję krawędzie wychodzące, przy czym bloki są już odpowiednio 
zorientowane. Metodzie add_out_edges nie dawałam za argument typu krawędzi, bo skoro dana krawędź 
wychodzi z danego bloku, a bloki, które ta krawedź łączy są już dobrze zorientowane, to typ krawędzi 
jest zawsze (1, -1). Jednak teraz gubię informacje, których krawędzi nie wprowadziłam do grafu :( Można
by zwrócić listę takich krawędzi. Masz może jakis lepszy pomysł jak rozwiązać ten problem? """            
        

        
            
        
        