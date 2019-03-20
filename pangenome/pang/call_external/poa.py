import os
from subprocess import run
from pathlib import Path

from tools import loggingtools

detailed_logger = loggingtools.get_logger("details")
def call(po_file_path: Path, hb_file_path: Path, blosum_path: Path, hbmin:float) -> None:
    poa_path = Path(os.path.abspath(__file__)).joinpath('../../bin/poa').resolve()
    detailed_logger.info(f"Run poa! Input: {po_file_path} Output: {hb_file_path}...")
    run([poa_path, '-read_msa', po_file_path, '-hb', '-po', hb_file_path, blosum_path,  '-hbmin',
         str(hbmin)])
