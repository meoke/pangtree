from subprocess import run

def remove_shifted_or_reversed(mycoplasma_file_name):
    run(['sed', '\'/GCF_000733995\|GCF_000387745/d\|label=u0\|label=u1\|label=u2\'', mycoplasma_file_name])
