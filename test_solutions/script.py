import timeit
from collections import namedtuple
from SequenceInfo import SequenceInfo
from EdgeType import EdgeType

# Możliwe rozwiązania:

# 1. słownik, w którym klucz to namedtuple
Edge = namedtuple('Edge', ['left', 'right', 'edgetype'])
edges1 = []
for i in range(10):
    edges1.append(Edge(i, i+1, EdgeType.LEFT_TO_LEFT))
    edges1.append(Edge(i, i+1, EdgeType.LEFT_TO_RIGHT))
    edges1.append(Edge(i, i+1, EdgeType.RIGHT_TO_LEFT))
    edges1.append(Edge(i, i+1, EdgeType.LEFT_TO_RIGHT))

edges_dict1 = {edge: ([], []) for edge in edges1}

# 1. słownik, w którym klucz to namedtuple,  a edge type to zwykły int
edges1prim = []
for i in range(10):
    edges1prim.append(Edge(i, i+1, 0))
    edges1prim.append(Edge(i, i+1, 1))
    edges1prim.append(Edge(i, i+1, 2))
    edges1prim.append(Edge(i, i+1, 3))

edges_dict1prim = {edge: ([], []) for edge in edges1prim}

# 2. słownik, w którym klucz to tuple
edges2 = []
for i in range(10):
    edges2.append((i, i + 1, EdgeType.LEFT_TO_LEFT))
    edges2.append((i, i + 1, EdgeType.LEFT_TO_RIGHT))
    edges2.append((i, i + 1, EdgeType.RIGHT_TO_LEFT))
    edges2.append((i, i + 1, EdgeType.LEFT_TO_RIGHT))

edges_dict2 = {edge: ([], []) for edge in edges2}

# 2prim. słownik, w którym klucz to tuple, a edge type to zwykły int
edges2prim = []
for i in range(10):
    edges2prim.append((i, i + 1, 0))
    edges2prim.append((i, i + 1, 1))
    edges2prim.append((i, i + 1, 2))
    edges2prim.append((i, i + 1, 3))

edges_dict2prim = {edge: ([], []) for edge in edges2prim}

# 3. słownik słowników słowników :D
edges3 = {}
for i in range(10):
    edges3[i] = {}
    for j in range(10):
        edges3[i][j] = {EdgeType.RIGHT_TO_LEFT: ([], []),
                        EdgeType.LEFT_TO_RIGHT: ([], []),
                        EdgeType.RIGHT_TO_LEFT: ([], []),
                        EdgeType.RIGHT_TO_RIGHT: ([], [])}


def solution1():
    e = Edge(8, 9, EdgeType.RIGHT_TO_LEFT)
    edges_dict1[e][0].append(SequenceInfo(0, 0))


def solution1prim():
    e = Edge(8, 9, 0)
    edges_dict1prim[e][0].append(SequenceInfo(0, 0))


def solution2():
    e = (8, 9, EdgeType.RIGHT_TO_LEFT)
    edges_dict2[e][0].append(SequenceInfo(0, 0))


def solution2prim():
    e = (8, 9, 0)
    edges_dict2prim[e][0].append(SequenceInfo(0, 0))


def solution3():
    edges3[8][9][EdgeType.RIGHT_TO_LEFT][0].append(SequenceInfo(0, 0))


print(f"Solution 1: {timeit.timeit('solution1()', globals=globals())}")
print(f"Solution 1prim: {timeit.timeit('solution1prim()', globals=globals())}")
print(f"Solution 2: {timeit.timeit('solution2()', globals=globals())}")
print(f"Solution 2prim: {timeit.timeit('solution2prim()', globals=globals())}")
print(f"Solution 3: {timeit.timeit('solution3()', globals=globals())}")


#uwaga, wykorzystanie po prostu 0, 1, 2, 3 zamiast tego enuma EdgeTypa powoduje prawie dwukrotne przyspieszenie dla solution2 (i w innych też jest znaczące przyspieszenie