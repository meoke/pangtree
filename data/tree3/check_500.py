from Bio import SeqIO
to_save = []
ids = ["AB", "CD", "ABCD", "EFGHIJ", "GHIJ", "EF", "HIJ"]
with open("all_500_corrected.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if record.id in ids:
            to_save.append((f"{record.id}", record.seq))

for p in to_save:
    with open("fa_internal_nodes/" + p[0] + ".fasta", "w") as output:
        output.writelines(f">{p[0]}\n")
        output.writelines(f"{p[1]}\n")