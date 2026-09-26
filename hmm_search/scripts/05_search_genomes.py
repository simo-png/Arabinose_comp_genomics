#!/usr/bin/env python3
"""Runs hmmsearch for the combined gene-profile database (built by 04_press_db.sh)
against every genome .faa file in a user-supplied directory. One hmmsearch
invocation per genome: hmmsearch iterates all models in the combined profile
file internally, so this is 1 run per genome, not 1 run per gene per genome.
Output: results/raw/<genome_stem>.tblout, plus a matching log per genome.
"""

import argparse
import subprocess
import sys
from pathlib import Path

CPU = 10
OUT_DIR = "Streptomycetae"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "genome_faa_dir",
        type=Path,
        help="Directory containing one .faa protein file per genome",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    genome_dir: Path = args.genome_faa_dir

    if not genome_dir.is_dir():
        print(f"Genome directory not found: {genome_dir}", file=sys.stderr)
        return 1

    script_dir = Path(__file__).resolve().parent
    project_dir = script_dir.parent

    hmm_db = project_dir / "profiles" / "streptomycetae_conserved.hmm"
    out_dir = project_dir / "results" / "raw" / OUT_DIR
    log_dir = project_dir / "logs"

    if not hmm_db.is_file():
        print(f"Combined profile database not found: {hmm_db}", file=sys.stderr)
        print("Run 04_press_db.sh first.", file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    faa_files = sorted(genome_dir.glob("*.faa"))
    if not faa_files:
        print(f"No .faa files found in {genome_dir}", file=sys.stderr)
        return 1

    print(f"Searching {len(faa_files)} genomes in {genome_dir} against {hmm_db.name}")

    for faa_file in faa_files:
        stem = faa_file.stem
        out_file = out_dir / f"{stem}.tblout"
        log_file = log_dir / f"{stem}_hmmsearch.log"

        print(f"Searching {faa_file.name} -> {out_file.name}")

        with open(log_file, "w") as log_fh:
            subprocess.run(
                [
                    "hmmsearch",
                    "--cpu",
                    str(CPU),
                    "--tblout",
                    str(out_file),
                    str(hmm_db),
                    str(faa_file),
                ],
                stdout=log_fh,
                stderr=subprocess.STDOUT,
                check=True,
            )

    print(f"Done. Per-genome tblout results written to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
