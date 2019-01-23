from dashsite.pang.MafgraphAnia import sorter

maf_path = "../examples/Mycoplasma/alignment.maf"
maf_path_ania = "../examples/blocks/mulAli.maf"
maf_path_paulina = "../examples/blocks/smallMulAli.maf"
maf_path_paulina_b1 = "../examples/blocks/b1.maf"
maf_path_paulina_b1b = "../examples/blocks/b1B.maf"
maf_path_paulina_b2 = "../examples/blocks/b2.maf"
blocks = sorter.sort_mafblocks(maf_path_paulina_b1b)
pass