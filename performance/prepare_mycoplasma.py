from Bio import AlignIO

myco_path = '../examples/m/alignment_clean.txt'
myco_1 = "../examples/m/first_block.maf"
seqrecs = set()
# p = [*AlignIO.parse(myco_path, "maf")][0]
# AlignIO.write(p, myco_1, "maf" )
for multiple_alignment in AlignIO.parse(myco_1, "maf"):
    print("printing a new multiple alignment")

    for seqrec in multiple_alignment:
        seqrecs.add(seqrec.id)
        print("%s starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
              (seqrec.id,
               seqrec.annotations["start"],
               seqrec.annotations["strand"],
               seqrec.annotations["srcSize"],
               seqrec.annotations["size"]))

# for multiple_alignment in AlignIO.parse(myco_path, "maf"):
#     multiple_alignment.
#     print("printing a new multiple alignment")
#
#     for seqrec in multiple_alignment:
#         seqrecs.add(seqrec.id)
#         # print("%s starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
#         #       (seqrec.id,
#         #        seqrec.annotations["start"],
#         #        seqrec.annotations["strand"],
#         #        seqrec.annotations["srcSize"],
#         #        seqrec.annotations["size"]))
# print("\n".join(list(seqrecs)))

# with open(myco_path) as myco:
#     mycolines = myco.readlines()
#     print(len(mycolines))