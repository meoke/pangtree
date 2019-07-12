
name = "tree3/tree3_100_nocycles"
ext = ".maf"
# with open(name+ext) as input:
#     i = input.read()
#     i = i.replace("A.chr20", "chr20.A")
#     i = i.replace("B.chr20", "chr20.B")
#     i = i.replace("C.chr20", "chr20.C")
#     i = i.replace("D.chr20", "chr20.D")
#     i = i.replace("E.chr20", "chr20.E")
#     i = i.replace("F.chr20", "chr20.F")
#     i = i.replace("G.chr20", "chr20.G")
#     i = i.replace("H.chr20", "chr20.H")
#     i = i.replace("I.chr20", "chr20.I")
#     i = i.replace("J.chr20", "chr20.J")
#
# with open(name+"_corrected"+ext, "w") as output:
#     output.write(i)

with open(name + ext) as input:
    i = input.readlines()
    for l in i:
        if "I.chr20" in l:
            print(l)


# from Bio import SeqIO
# with open("example.fasta", "rU") as handle:
#     for record in SeqIO.parse(handle, "fasta"):
#         print(record.id)