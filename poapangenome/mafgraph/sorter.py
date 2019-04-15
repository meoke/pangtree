# -*- coding: utf-8 -*-
import networkx as nx
import six
from .mafreader import read_maf
from .weighted_edges import weight, Edge

def _out_edges(v, blocks, ub, G):
    # Return a list of nodes connected to node v by edge leaving the vertex v   
    return [u for u in G[v] if ub >= blocks[u].order() > blocks[v].order()]

def _in_edges(v, blocks, lb, G):
    # Return a list of nodes connected to node v by edge coming in to the vertex v     
    return [u for u in G[v] if lb < blocks[u].order() < blocks[v].order()]  

def _dfs_f(root, ub, G, blocks):
    visited = set()
    stack = [root,]
    while stack:
        node = stack.pop()
        if blocks[node].order()==ub: return []
        if node not in visited:
            visited.add(node)
            stack.extend([x for x in _out_edges(node, blocks, ub, G) if x not in visited])
    return list(visited) 

def _dfs_b(root, lb, G, blocks):
    visited = set()
    stack = [root,]
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            stack.extend([x for x in _in_edges(node, blocks, lb, G) if x not in visited])
    return list(visited)

def _reorder(R_f, R_b, blocks):
    R_f.sort(key = lambda x: blocks[x].order()) 
    R_b.sort(key = lambda x: blocks[x].order())
    L = R_b + R_f
    O = sorted(blocks[x].order() for x in L)
    if six.PY2:
        for i in xrange(len(L)):
            blocks[L[i]].reorder(O[i])
    else:
        for i in range(len(L)):
            blocks[L[i]].reorder(O[i])

def _add_edge_within_component(x, y, G, blocks):
    # Add edge between blocks of the same component
    lb = blocks[y].order()
    ub = blocks[x].order()
    if lb is ub: return
    elif lb < ub:
        R_f = _dfs_f(y, ub, G, blocks)
        if R_f:
            R_b = _dfs_b(x, lb, G, blocks)
            _reorder(R_f, R_b, blocks)
            G.add_edge(x,y)
    else: 
        G.add_edge(x,y)
       
def _add_edge_between_components(e, blocks):
    # Add edge between blocks of different components
    reverse, flank = 1, 1
    if blocks[e.left].size() < blocks[e.right].size():
        e = Edge(e.right, e.left, e.type[::-1])
    if blocks[e.left].orientation()*blocks[e.right].orientation() is e.type[0]*e.type[1]:
        reverse = -1
    if blocks[e.left].orientation()*e.type[0] < 0: 
        flank = -1
    blocks[e.right].unionto(blocks[e.left], reverse, flank)

def connect_components(blocks):
    if len(blocks) != blocks[0].size(): 
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
            block.orient_block()
    else:
        for block in blocks:
            block.orient_block()
            
def set_out_edges(d, blocks):
    for edge in d:
        edge_type = (blocks[edge.left].orientation()*edge.type[0], blocks[edge.right].orientation()*edge.type[1])
        for tup in d[edge][0]:
            tup[0].set_start_position(blocks[edge.left].alignment, blocks[edge.left].orientation())
            tup[1].set_start_position(blocks[edge.right].alignment, blocks[edge.right].orientation())   
        if blocks[edge.left].order() < blocks[edge.right].order():
            blocks[edge.left].add_out_edges(edge.right, edge_type, d[edge][0])
        elif blocks[edge.left].order() > blocks[edge.right].order():
            sequences = [x[::-1] for x in d[edge][0]]
            edge_type = edge_type[::-1]
            blocks[edge.right].add_out_edges(edge.left, edge_type, sequences)
        else:
            blocks[edge.left].add_out_edges(edge.right, edge.type, d[edge][0])

def sort_mafblocks(maf_file):   
    blocks, seq = read_maf(maf_file) # blocks - list of Block instances
    d = weight(seq)
    edges = sorted(d.keys(), key=lambda x: (d[x][1], x.type, x.left, x.right)) # list of edges sorted by the weight
    G = nx.Graph()
    for e in edges:
        if blocks[e.left].find() is blocks[e.right].find():
            if blocks[e.left].orientation()*blocks[e.right].orientation() is e.type[0]*e.type[1]:
                continue
            elif blocks[e.left].orientation()*e.type[0] > 0:
                _add_edge_within_component(e.left, e.right, G, blocks)
            else:
                _add_edge_within_component(e.right, e.left, G, blocks)           
        else:
            G.add_edge(e.left, e.right)
            _add_edge_between_components(e, blocks)
    set_out_edges(d, blocks)
    connect_components(blocks)
    blocks = sorted(blocks, key = lambda b: b.order())
    return blocks

