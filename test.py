# from mafgraph import sorter
#
# b = sorter.sort_mafblocks('examples/blocks/b2.maf')
# pass 739080987

from Bio import AlignIO
ebola_path = 'examples/Ebola/ebola.maf'
path = "/home/meoke/PycharmProjects/pang/tests/PangraphBuilder_Tests/PangraphBuilderFromMAF_Tests/files_build_pangraph/test_1_messy_sequences.maf"
a = AlignIO.parse(path, "maf")
pass
# size = 0
# d = {}
# for multiple_alignment in AlignIO.parse(ebola_path, "maf"):
#     for seqrec in multiple_alignment:
#         if seqrec.id not in d.keys():
#             d[seqrec.id] = (seqrec.annotations["size"], seqrec.annotations["srcSize"])
#         else:
#             d[seqrec.id] = (d[seqrec.id][0] + seqrec.annotations["size"], d[seqrec.id][1])
# for k, v in d.items():
#     if v[0] != v[1]:
#         print(f"Seq {k}, \t\texpected: {v[1]}, \t\tactual: {v[0]}")

# from Bio import Entrez
# Entrez.email = "pedziadkiewicz@gmail.com"  # Always tell NCBI who you are
# # handle = Entrez.esearch(db="nucleotide", accession="KC242801", version="KC242801.1")
# handle = Entrez.efetch(db="nucleotide", id="KJ660346.2", rettype="fasta", retmode="text", seq_start=0, seq_stop=10)
# record = handle.read()
# print(record)
# handle.close()
id1 = "NC_024781v1.NC_024781v1"
id2 = "NC_001608v3.NC_001608v3"
# ??? To: 47 Node_id: 1611, Seq_id: NC_024781v1.NC_024781v1, Seq pos: 153
# ??? To: 47 Node_id: 1611, Seq_id: NC_001608v3.NC_001608v3, Seq pos: 153

for i, multiple_alignment in enumerate(a):
    for seqrec in multiple_alignment:
        if seqrec.id == id2:
            print("%s starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
              (i, seqrec.annotations["start"],
               seqrec.annotations["strand"],
               seqrec.annotations["srcSize"],
               seqrec.annotations["size"]))