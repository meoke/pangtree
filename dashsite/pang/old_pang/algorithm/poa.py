from pathlib import Path
from call_external import poa_runner


def run(input_path: Path, output_path: Path, hbmin: float) -> None:
    poa_runner.run_poa(po_file_path=input_path,
                       hb_file_path=output_path,
                       hbmin=hbmin)
