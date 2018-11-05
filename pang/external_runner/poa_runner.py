from subprocess import run
from pathlib import Path

def run_poa(po_file_path: Path, hb_file_path: Path, hbmin):
    poa_path = po_file_path.joinpath('../../../bin/poa').resolve()
    blosum_path = po_file_path.joinpath('../../../bin/blosum80.mat').resolve()
    run([poa_path, '-read_msa', po_file_path, '-hb', '-po', hb_file_path, blosum_path,  '-hbmin',
         str(hbmin)])