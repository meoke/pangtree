import os
import subprocess
from pathlib import Path

from pangenome.pang.tools import loggingtools

detailed_logger = loggingtools.get_logger("details")


def call(po_file_path: Path, hb_file_path: Path, blosum_path: Path, hbmin: float) -> None:
    poa_path = Path(os.path.abspath(__file__)).joinpath('../../bin/poa').resolve()
    detailed_logger.info(f"Run poa! Input: {po_file_path} Output: {hb_file_path}...")
    command = f"{poa_path} -read_msa {po_file_path} -hb -po {hb_file_path} {blosum_path} -hbmin {hbmin}"
    poa_result = subprocess.run(command, stderr=subprocess.PIPE, shell=True)
    poa_str_output = poa_result.stderr.decode("ASCII")
    detailed_logger.info(f"Poa output: {poa_str_output}")
