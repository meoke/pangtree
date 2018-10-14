import sys
from Bio import AlignIO
from maf2.Mafgraph.Mafgraph import Mafgraph


mafpath = sys.argv[1]
maf = AlignIO.parse(mafpath, "maf")
mafgraph = Mafgraph(maf, remove_cycles=False)
print(mafgraph)
