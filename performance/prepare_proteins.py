from Bio import AlignIO

myco_path = '../examples/proteins/proteins.maf'
seqrecs = set()

def calc_letters(s):
    c = 0
    for i in s:
        if i != '-':
            c += 1
    return c

for multiple_alignment in AlignIO.parse(myco_path, "maf"):
    for seqrec in multiple_alignment:
        l = calc_letters(seqrec.seq)#seqrec.annotations["start"])
        print(f"{seqrec.id}: {l}")
        # print("%s starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
        #       (seqrec.id,
        #        seqrec.annotations["start"],
        #        seqrec.annotations["strand"],
        #        seqrec.annotations["srcSize"],
        #        seqrec.annotations["size"]))



print(calc_letters("---AA"))
print(calc_letters("---"))
print(calc_letters("AAA"))