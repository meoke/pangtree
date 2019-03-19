import os
from subprocess import run
from pathlib import Path
import logging


def call(po_file_path: Path, hb_file_path: Path, blosum_path: Path, hbmin:float) -> None:
    poa_path = Path(os.path.abspath(__file__)).joinpath('../../bin/poa').resolve()
    logging.info(f"Run poa on {po_file_path} and produce {hb_file_path}...")
    run([poa_path, '-read_msa', po_file_path, '-hb', '-po', hb_file_path, blosum_path,  '-hbmin',
         str(hbmin)])
