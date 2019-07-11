
name = "all_500"
ext = ".fasta"
with open(name+ext) as input:
    i = input.read()
    i = i.replace(".chr20", "")

with open(name+"_corrected"+ext, "w") as output:
    output.write(i)